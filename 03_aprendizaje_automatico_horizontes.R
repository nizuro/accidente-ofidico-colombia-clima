# =============================================================================
# Accidente ofídico en Colombia, 2015-2024
# Script 03. Aprendizaje automático y horizontes de predicción (1, 4 y 8 semanas)
#
# Requiere: 00_configuracion.R
# Entrada : dataset_modelado_brotes.rds (salida del script 01)
# Salida  : resultados_R/tabla_conteo_horizontes.csv, tabla_brote_horizontes.csv,
#           importancia_conteo.csv, importancia_rf_brote.csv y las figuras g07 a g09
#
# Tarea A (conteo, municipio-semana). Modelo de Poisson con intercepto por
#   municipio y estacionalidad propia de cada departamento (estadístico) frente
#   a boosting de Poisson (aprendizaje automático), con y sin variables climáticas.
# Tarea B (brote, departamento-semana). Regresión logística frente a bosque
#   aleatorio y boosting, con y sin variables climáticas.
#
# Regla anti-fuga: para predecir la semana t con horizonte h solo se usa
# información hasta la semana t-h. La estacionalidad (semana del año de t) y la
# población se consideran conocidas.
#
# Notas de ejecución
#   - Cada horizonte guarda un punto de control (ckpt_*.rds). Si el proceso se
#     interrumpe, al volver a ejecutarlo se retoma donde quedó. Para recalcular
#     todo, hay que borrar esos archivos o poner REUSAR_CHECKPOINTS <- FALSE.
#   - El boosting usa el paquete gbm y el bosque aleatorio el paquete ranger.
#   - Los gráficos se dibujan en el panel Plots de RStudio y se guardan en PNG.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Configuración
# -----------------------------------------------------------------------------
nucleo <- c("dplyr", "readr", "tidyr", "ggplot2", "mgcv")
faltan <- nucleo[!sapply(nucleo, requireNamespace, quietly = TRUE)]
if (length(faltan)) install.packages(faltan)
invisible(lapply(nucleo, library, character.only = TRUE))

source("00_configuracion.R")

# gbm y ranger se instalan, si faltan, en su versión binaria en Windows (evita compilar).
for (p in c("gbm", "ranger"))
  if (!instalar_bin(p)) stop("No se pudo instalar el paquete ", p, ": install.packages('", p, "')")

AJUSTAR_BAM        <- TRUE         # FALSE omite los modelos aditivos estadísticos (mucho más rápido)

# -----------------------------------------------------------------------------
# 1. Datos y funciones comunes
# -----------------------------------------------------------------------------
df <- readRDS(file.path(DIR_DATOS, "dataset_modelado_brotes.rds"))

# Anomalías climáticas contemporáneas, rezagadas una semana (disponibles hasta t-1)
vars_contemp <- c("z_precip_4w", "z_precip_13w", "z_precip_26w", "z_temp_4w", "z_rh_4w")
df <- df |>
  arrange(municipio_id, week_start) |>
  group_by(municipio_id) |>
  mutate(across(all_of(vars_contemp), ~ dplyr::lag(.x, 1), .names = "{.col}_l1")) |>
  ungroup() |>
  mutate(features_completas = as.logical(features_completas),
         pob_rural_prop = poblacion_rural / poblacion_total)

# Columnas disponibles hasta t-1. Para un horizonte h se desplazan h-1 semanas más.
auto_cols  <- c("casos_lag1", "casos_lag2", "casos_lag4", "casos_sum4w_prev",
                "casos_sum13w_prev", "casos_dpto_sum4w_prev")
clima_cols <- c("z_precip_4w_l1", "z_precip_13w_l1", "z_precip_26w_l1",
                "z_temp_4w_l1", "z_rh_4w_l1",
                "precip_mm_lag1", "precip_mm_lag2", "precip_mm_lag4", "precip_mm_lag8",
                "temp_mean_lag1", "rh_mean_lag1")
fijas_cols <- c("log_poblacion_total", "pob_rural_prop", "sin_semana", "cos_semana")

resumen <- function(a, nivel) {
  tibble(nivel = nivel, MAE = mean(abs(a$obs - a$pred)),
         RMSE = sqrt(mean((a$obs - a$pred)^2)),
         R2 = 1 - sum((a$obs - a$pred)^2) / sum((a$obs - mean(a$obs))^2),
         dev_poisson = pois_dev(a$obs, a$pred),
         casos_obs = sum(a$obs), casos_pred = round(sum(a$pred)))
}

# Métricas de conteo en tres niveles de agregación: municipio, departamento y nación
eval_conteo <- function(te, pred, modelo, h) {
  d <- te |> select(municipio_id, departamento, week_start, casos) |> mutate(pred = pred)
  bind_rows(
    resumen(d |> group_by(municipio_id, week_start) |>
              summarise(obs = sum(casos), pred = sum(pred), .groups = "drop"), "municipio-semana"),
    resumen(d |> group_by(departamento, week_start) |>
              summarise(obs = sum(casos), pred = sum(pred), .groups = "drop"), "departamento-semana"),
    resumen(d |> group_by(week_start) |>
              summarise(obs = sum(casos), pred = sum(pred), .groups = "drop"), "nacional-semana")
  ) |> mutate(modelo = modelo, horizonte = h, .before = 1)
}


# -----------------------------------------------------------------------------
# 2. Tarea A. Conteo municipal por horizonte
# -----------------------------------------------------------------------------
# Datos de entrenamiento y prueba para un horizonte h. Las variables que se
# conocen hasta t-1 se desplazan h-1 semanas más, para que solo se use información
# hasta t-h. Los municipios sin historia de entrenamiento se retiran de la prueba.
prep_mun <- function(h) {
  d <- df
  if (h > 1) {
    d <- d |> group_by(municipio_id) |>
      mutate(across(all_of(c(auto_cols, clima_cols)), ~ dplyr::lag(.x, h - 1))) |>
      ungroup()
  }
  d <- d |>
    filter(features_completas, poblacion_total > 0, !is.na(log_poblacion_total)) |>
    drop_na(all_of(c(auto_cols, clima_cols, fijas_cols))) |>
    mutate(departamento = factor(departamento), municipio_id = factor(municipio_id))
  tr <- d |> filter(particion == "entrenamiento") |> droplevels()
  te <- d |> filter(particion == "prueba", municipio_id %in% levels(tr$municipio_id))
  te$municipio_id <- factor(te$municipio_id, levels = levels(tr$municipio_id))
  te$departamento <- factor(te$departamento, levels = levels(tr$departamento))
  list(tr = tr, te = te)
}

# Estacionalidad propia de cada departamento: mejoró el AIC en unos 280 puntos
# frente a una curva común nacional (bloque 3e del script 02).
f_base <- casos ~ offset(log_poblacion_total) + departamento +
  s(semana_anio, by = departamento, bs = "cc", k = 8) +
  log1p(casos_lag1) + log1p(casos_sum4w_prev) + log1p(casos_sum13w_prev) +
  log1p(casos_dpto_sum4w_prev) + s(municipio_id, bs = "re")

# Modelo con clima: precipitación (anomalías y lluvia semanal rezagada), temperatura y humedad
f_clim <- update(f_base, . ~ . + s(z_precip_4w_l1, k = 6) + s(z_precip_13w_l1, k = 6) +
                   s(z_precip_26w_l1, k = 6) + log1p(precip_mm_lag1) + log1p(precip_mm_lag2) +
                   log1p(precip_mm_lag4) + log1p(precip_mm_lag8) +
                   s(z_temp_4w_l1, k = 6) + s(z_rh_4w_l1, k = 6))

# Boosting de Poisson con offset de población. Devuelve la predicción y, en el
# modelo con clima, la influencia relativa de cada variable.
ml_conteo <- function(tr, te, usar_clima) {
  vars <- c(auto_cols, fijas_cols, if (usar_clima) clima_cols)
  f <- as.formula(paste("casos ~ offset(log_poblacion_total) +",
                        paste(c(vars, "departamento"), collapse = " + ")))
  set.seed(SEMILLA)   # gbm usa bag.fraction: fijar la semilla hace reproducible el ajuste
  m <- gbm::gbm(f, data = as.data.frame(tr), distribution = "poisson", n.trees = 300,
                interaction.depth = 4, shrinkage = 0.05, bag.fraction = 0.8,
                n.minobsinnode = 50, verbose = FALSE)
  # predict.gbm no suma el offset: se agrega a mano en la escala logarítmica.
  pred <- exp(predict(m, as.data.frame(te), n.trees = 300, type = "link") + te$log_poblacion_total)
  imp <- NULL
  if (usar_clima) {
    s <- summary(m, plotit = FALSE)
    imp <- setNames(s$rel.inf, as.character(s$var))
  }
  list(pred = as.numeric(pred), imp = imp, metodo = "influencia relativa (%)")
}

res_conteo <- list(); imp_conteo <- list()
for (h in HORIZONTES) {
  ckpt <- file.path(DIR_SALIDA, sprintf("ckpt_A_h%d.rds", h))
  if (REUSAR_CHECKPOINTS && file.exists(ckpt)) {
    cat(sprintf("\nTarea A, horizonte %d: recuperado del punto de control\n", h))
    r <- readRDS(ckpt)
    res_conteo[[length(res_conteo) + 1]] <- r$res
    imp_conteo[[length(imp_conteo) + 1]] <- r$imp
    next
  }
  cat(sprintf("\nTarea A, horizonte de %d semana(s)\n", h))
  p <- prep_mun(h); tr <- p$tr; te <- p$te
  cat(sprintf("Entrenamiento: %d filas | Prueba: %d filas\n", nrow(tr), nrow(te)))
  res_h <- list()

  # Referencia: tasa histórica del municipio por su población
  tasa_m <- tr |> group_by(municipio_id) |>
    summarise(tasa = sum(casos) / sum(poblacion_total), .groups = "drop")
  p_media <- te |> left_join(tasa_m, by = "municipio_id") |>
    mutate(p = coalesce(tasa, 0) * poblacion_total) |> pull(p)
  res_h[["media"]] <- eval_conteo(te, p_media, "media_municipal", h)

  if (AJUSTAR_BAM) {
    cat("Poisson con intercepto por municipio y estacionalidad por departamento...\n")
    aic_auto <- NA
    for (nm in c("pois_re_auto", "pois_re_clima")) {
      f <- if (nm == "pois_re_auto") f_base else f_clim
      m <- mgcv::bam(f, family = poisson(), data = tr, method = "fREML",
                     discrete = TRUE, nthreads = NTHREADS)
      pr <- as.numeric(predict(m, newdata = te, type = "response", n.threads = NTHREADS))
      res_h[[nm]] <- eval_conteo(te, pr, nm, h)
      if (nm == "pois_re_auto") aic_auto <- AIC(m) else
        cat(sprintf("  AIC sin clima %.1f | con clima %.1f | diferencia %.1f\n",
                    aic_auto, AIC(m), AIC(m) - aic_auto))
    }
  }

  cat("Boosting de Poisson, sin y con clima...\n")
  imp_h <- NULL
  for (usar in c(FALSE, TRUE)) {
    fit <- ml_conteo(tr, te, usar)
    nm <- if (usar) "boost_clima" else "boost_auto"
    res_h[[nm]] <- eval_conteo(te, fit$pred, nm, h)
    if (usar && !is.null(fit$imp))
      imp_h <- tibble(horizonte = h, variable = names(fit$imp),
                      importancia = as.numeric(fit$imp), metodo = fit$metodo)
  }
  res_h <- bind_rows(res_h)
  saveRDS(list(res = res_h, imp = imp_h), ckpt)
  res_conteo[[length(res_conteo) + 1]] <- res_h
  imp_conteo[[length(imp_conteo) + 1]] <- imp_h
}

tabla_conteo <- bind_rows(res_conteo) |> mutate(across(where(is.numeric), ~ round(.x, 3)))
cat("\nConteo: métricas en el periodo de prueba\n")
print(tabla_conteo |> filter(nivel != "nacional-semana") |>
        arrange(horizonte, nivel, dev_poisson), n = 100)
write_csv(tabla_conteo, file.path(DIR_SALIDA, "tabla_conteo_horizontes.csv"))

# Desviación de Poisson por horizonte y modelo (menor es mejor)
g9 <- tabla_conteo |> filter(nivel != "nacional-semana") |>
  ggplot(aes(factor(horizonte), dev_poisson, color = modelo, group = modelo)) +
  geom_line() + geom_point(size = 2) + facet_wrap(~ nivel, scales = "free_y") +
  labs(title = "Error de predicción de conteos según la anticipación",
       subtitle = "Desviación de Poisson en 2022-2024 (menor = mejor)",
       x = "Semanas de anticipación", y = "Desviación de Poisson", color = NULL) + tema
mostrar(g9, "g09_error_conteo_por_horizonte.png", 10, 4.5)

imp_tab <- bind_rows(imp_conteo)
if (nrow(imp_tab)) {
  imp_tab <- imp_tab |> filter(!grepl("^departamento", variable)) |>
    arrange(horizonte, desc(importancia))
  write_csv(imp_tab, file.path(DIR_SALIDA, "importancia_conteo.csv"))
  cat("\nVariables más importantes (boosting con clima):\n")
  print(imp_tab |> group_by(horizonte) |> slice_max(importancia, n = 6), n = 30)
  g8 <- imp_tab |> group_by(horizonte) |> slice_max(importancia, n = 10) |> ungroup() |>
    mutate(var_h = paste(variable, horizonte, sep = "___")) |>
    ggplot(aes(importancia, reorder(var_h, importancia))) +
    geom_col(fill = "#2b6a99") +
    scale_y_discrete(labels = function(x) sub("___.*", "", x)) +
    facet_wrap(~ paste0("Horizonte ", horizonte, " sem."), scales = "free") +
    labs(title = "Importancia de variables (boosting con clima)",
         subtitle = paste("Medida:", imp_tab$metodo[1]), x = NULL, y = NULL) + tema
  mostrar(g8, "g08_importancia_variables.png", 12, 5)
}


# -----------------------------------------------------------------------------
# 3. Tarea B. Brote departamental (clasificación)
# -----------------------------------------------------------------------------
dep <- df |>
  group_by(departamento, week_start, anio, semana_anio, particion) |>
  summarise(casos = first(casos_dpto), umbral = first(umbral_brote_dpto),
            brote = first(brote_dpto),
            across(c(z_precip_4w, z_precip_13w, z_temp_4w, z_rh_4w), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

# Suma de x entre los rezagos a y b (a <= b)
rollprev <- function(x, a, b) Reduce(`+`, lapply(a:b, function(k) dplyr::lag(x, k)))

# Variables de la semana t construidas solo con información hasta t-h
prep_dep <- function(h) {
  dep |> arrange(departamento, week_start) |> group_by(departamento) |>
    mutate(c_lag = dplyr::lag(casos, h), c_sum4 = rollprev(casos, h, h + 3),
           c_sum13 = rollprev(casos, h, h + 12), brote_l = dplyr::lag(brote, h),
           zp4_l = dplyr::lag(z_precip_4w, h), zp13_l = dplyr::lag(z_precip_13w, h),
           zt4_l = dplyr::lag(z_temp_4w, h), zr4_l = dplyr::lag(z_rh_4w, h)) |>
    ungroup() |>
    drop_na(c_lag, c_sum4, c_sum13, brote_l, zp4_l, zp13_l, zt4_l, zr4_l, umbral, brote) |>
    mutate(departamento = factor(departamento),
           sin_s = sin(2 * pi * semana_anio / 52), cos_s = cos(2 * pi * semana_anio / 52),
           brote_f = factor(brote, levels = c(0, 1)))
}

# Área bajo la curva precisión-exhaustividad (precisión promedio)
pr_auc <- function(y, s) {
  o <- order(s, decreasing = TRUE); y <- y[o]
  sum((cumsum(y) / seq_along(y))[y == 1]) / sum(y)
}

# Métricas de clasificación. La alerta se emite para el 25 % de las semanas con
# mayor puntaje; la precisión y la exhaustividad se calculan sobre esa alerta.
metricas_brote <- function(y, s, modelo, h) {
  corte <- quantile(s, 0.75); alerta <- s >= corte
  tibble(horizonte = h, modelo = modelo, AUC = auc_rank(y, s), PR_AUC = pr_auc(y, s),
         Brier = mean((s - y)^2), precision_top25 = sum(alerta & y == 1) / sum(alerta),
         recall_top25 = sum(alerta & y == 1) / sum(y), prevalencia = mean(y))
}

vars_aut <- c("departamento", "c_lag", "c_sum4", "c_sum13", "brote_l", "umbral", "sin_s", "cos_s")
vars_cli <- c(vars_aut, "zp4_l", "zp13_l", "zt4_l", "zr4_l")

# Boosting de Bernoulli; devuelve la probabilidad de brote en la muestra de prueba
boost_brote <- function(tr, te, vars) {
  f <- as.formula(paste("brote ~", paste(vars, collapse = " + ")))
  set.seed(SEMILLA)   # gbm usa bag.fraction: fijar la semilla hace reproducible el ajuste
  m <- gbm::gbm(f, data = as.data.frame(tr), distribution = "bernoulli", n.trees = 250,
                interaction.depth = 3, shrinkage = 0.05, bag.fraction = 0.8,
                n.minobsinnode = 10, verbose = FALSE)
  predict(m, as.data.frame(te), n.trees = 250, type = "response")
}

res_brote <- list(); imp_rf <- list()
for (h in HORIZONTES) {
  ckpt <- file.path(DIR_SALIDA, sprintf("ckpt_B_h%d.rds", h))
  if (REUSAR_CHECKPOINTS && file.exists(ckpt)) {
    cat(sprintf("\nTarea B, horizonte %d: recuperado del punto de control\n", h))
    r <- readRDS(ckpt)
    res_brote[[length(res_brote) + 1]] <- r$res
    imp_rf[[length(imp_rf) + 1]] <- r$imp
    next
  }
  cat(sprintf("\nTarea B, horizonte de %d semana(s)\n", h))
  d <- prep_dep(h)
  tr <- d |> filter(particion == "entrenamiento"); te <- d |> filter(particion == "prueba")
  te$departamento <- factor(te$departamento, levels = levels(tr$departamento))
  cat(sprintf("Entrenamiento: %d (brotes %.1f%%) | Prueba: %d (brotes %.1f%%)\n",
              nrow(tr), 100 * mean(tr$brote), nrow(te), 100 * mean(te$brote)))
  res_h <- list(); imp_h <- NULL

  for (con_clima in c(FALSE, TRUE)) {
    vars <- if (con_clima) vars_cli else vars_aut
    suf  <- if (con_clima) "_clima" else "_auto"

    # (1) Regresión logística
    m_lg <- suppressWarnings(glm(as.formula(paste("brote ~", paste(vars, collapse = " + "))),
                                 family = binomial(), data = tr))
    s_lg <- suppressWarnings(predict(m_lg, newdata = te, type = "response"))
    res_h[[paste0("lg", suf)]] <- metricas_brote(te$brote, s_lg, paste0("logistica", suf), h)

    # (2) Bosque aleatorio, con importancia por permutación
    m_rf <- ranger::ranger(as.formula(paste("brote_f ~", paste(vars, collapse = " + "))),
                           data = tr, probability = TRUE, num.trees = 400, min.node.size = 10,
                           importance = "permutation", seed = SEMILLA, num.threads = NTHREADS)
    s_rf <- predict(m_rf, data = te)$predictions[, "1"]
    res_h[[paste0("rf", suf)]] <- metricas_brote(te$brote, s_rf, paste0("random_forest", suf), h)
    if (con_clima)
      imp_h <- tibble(horizonte = h, variable = names(m_rf$variable.importance),
                      importancia = as.numeric(m_rf$variable.importance))

    # (3) Boosting
    res_h[[paste0("bo", suf)]] <- metricas_brote(te$brote, boost_brote(tr, te, vars),
                                                 paste0("boosting", suf), h)
  }

  # Referencia sin aprendizaje: persistencia (¿hubo brote hace h semanas?). El
  # segundo término desempata con los casos de esa semana.
  res_h[["pers"]] <- metricas_brote(te$brote, te$brote_l + 1e-6 * te$c_lag, "persistencia", h)

  res_h <- bind_rows(res_h)
  saveRDS(list(res = res_h, imp = imp_h), ckpt)
  res_brote[[length(res_brote) + 1]] <- res_h
  imp_rf[[length(imp_rf) + 1]] <- imp_h
}

tabla_brote <- bind_rows(res_brote) |> mutate(across(where(is.numeric), ~ round(.x, 3)))
cat("\nBrote: métricas en el periodo de prueba\n")
print(tabla_brote |> arrange(horizonte, desc(AUC)), n = 100)
write_csv(tabla_brote, file.path(DIR_SALIDA, "tabla_brote_horizontes.csv"))
if (length(imp_rf)) write_csv(bind_rows(imp_rf), file.path(DIR_SALIDA, "importancia_rf_brote.csv"))

g7 <- tabla_brote |>
  mutate(familia = sub("_(auto|clima)$", "", modelo),
         clima = if_else(grepl("_clima$", modelo), "con clima", "sin clima")) |>
  filter(familia != "persistencia") |>
  ggplot(aes(horizonte, AUC, color = familia, linetype = clima)) +
  geom_line() + geom_point() + scale_x_continuous(breaks = HORIZONTES) +
  geom_hline(yintercept = 0.5, color = "grey60") +
  labs(title = "Detección de brote departamental: AUC según anticipación",
       x = "Semanas de anticipación", y = "AUC (periodo de prueba 2022-2024)",
       color = NULL, linetype = NULL) + tema
mostrar(g7, "g07_auc_por_horizonte.png", 8, 4.5)

cat("\nScript 03 finalizado. Tablas y figuras en:", DIR_SALIDA, "\n")
