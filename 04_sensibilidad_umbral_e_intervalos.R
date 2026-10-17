# =============================================================================
# Accidente ofídico en Colombia, 2015-2024
# Script 04. Sensibilidad a la definición de brote e intervalos de confianza
#
# Requiere: 00_configuracion.R
# Entrada : dataset_modelado_brotes.rds (salida del script 01)
# Salida  : resultados_R/tabla_sensibilidad_umbral.csv, tabla_diferencias_auc.csv
#           y g10_sensibilidad_umbral.png
#
# Pregunta: ¿los resultados dependen del umbral con que se define un brote?
# Se repite la Tarea B (brote departamental a 1, 4 y 8 semanas) con tres
# definiciones, calculadas sobre el canal endémico de cada departamento:
#   p75      : casos > percentil 75   (definición usada en los demás scripts)
#   p90      : casos > percentil 90   (brotes más raros y más intensos)
#   media2de : casos > media + 2 desviaciones estándar
# En todas: años base 2015-2021 sin 2020, ventana de ±2 semanas y mínimo de 3 casos.
#
# Además se calculan:
#   - IC 95 % del AUC por bootstrap de bloques de semanas. Se remuestrean semanas
#     completas para conservar la correlación entre departamentos de una misma semana.
#   - IC 95 % de la diferencia de AUC entre modelos (con clima menos sin clima, y
#     modelo menos persistencia). Si el intervalo incluye el 0, no hay evidencia de
#     diferencia entre ambos.
#
# Cada combinación definición-horizonte guarda un punto de control
# (ckpt_S_<definicion>_h<h>.rds); para recalcular, hay que borrarlos.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Configuración
# -----------------------------------------------------------------------------
nucleo <- c("dplyr", "readr", "tidyr", "ggplot2")
faltan <- nucleo[!sapply(nucleo, requireNamespace, quietly = TRUE)]
if (length(faltan)) install.packages(faltan)
invisible(lapply(nucleo, library, character.only = TRUE))

source("00_configuracion.R")

# gbm y ranger se instalan, si faltan, en su versión binaria en Windows (evita compilar).
for (p in c("gbm", "ranger"))
  if (!instalar_bin(p)) stop("No se pudo instalar el paquete ", p, ": install.packages('", p, "')")

B_BOOT             <- 1000    # réplicas bootstrap

# -----------------------------------------------------------------------------
# 1. Tabla departamento-semana y etiquetas de brote
# -----------------------------------------------------------------------------
df <- readRDS(file.path(DIR_DATOS, "dataset_modelado_brotes.rds"))
dep <- df |>
  group_by(departamento, week_start, anio, semana_anio, particion) |>
  summarise(casos = first(casos_dpto),
            across(c(z_precip_4w, z_precip_13w, z_temp_4w, z_rh_4w), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

# Etiqueta de brote según la definición `def`. El estadístico del canal se calcula
# por departamento y semana del año, con ventana circular de ±VENTANA semanas.
etiquetar <- function(dep, def) {
  base <- dep |> filter(anio %in% ANIOS_BASE_CANAL, !anio %in% ANIOS_EXCLUIDOS_CANAL)
  canal <- base |> group_by(departamento) |>
    group_modify(~ {
      q <- sapply(1:53, function(s) {
        vec <- ((s - 1 + (-VENTANA_SEMANAS:VENTANA_SEMANAS)) %% 53) + 1
        x <- .x$casos[.x$semana_anio %in% vec]
        switch(def,
               p75      = unname(quantile(x, 0.75)),
               p90      = unname(quantile(x, 0.90)),
               media2de = mean(x) + 2 * sd(x))
      })
      tibble(semana_anio = 1:53, q = as.numeric(q))
    }) |> ungroup()
  dep |> left_join(canal, by = c("departamento", "semana_anio")) |>
    mutate(umbral = pmax(q, MIN_CASOS_BROTE - 1), brote = as.integer(casos > umbral))
}

cat("\nPrevalencia de brote por definición (% de semanas-departamento):\n")
prev <- bind_rows(lapply(DEFINICIONES, function(d) {
  e <- etiquetar(dep, d)
  tibble(definicion = d,
         prev_entrenamiento = mean(e$brote[e$particion == "entrenamiento"], na.rm = TRUE),
         prev_prueba = mean(e$brote[e$particion == "prueba"], na.rm = TRUE))
}))
print(prev |> mutate(across(where(is.numeric), ~ round(100 * .x, 1))))


# -----------------------------------------------------------------------------
# 2. Variables, modelos y métricas
# -----------------------------------------------------------------------------
# Suma de x entre los rezagos a y b (a <= b)
rollprev <- function(x, a, b) Reduce(`+`, lapply(a:b, function(k) dplyr::lag(x, k)))

# Variables de la semana t construidas solo con información hasta t-h
prep_dep <- function(depl, h) {
  depl |> arrange(departamento, week_start) |> group_by(departamento) |>
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

vars_aut <- c("departamento", "c_lag", "c_sum4", "c_sum13", "brote_l", "umbral", "sin_s", "cos_s")
vars_cli <- c(vars_aut, "zp4_l", "zp13_l", "zt4_l", "zr4_l")

# Boosting de Bernoulli; devuelve la probabilidad de brote en la muestra de prueba
boost_brote <- function(tr, te, vars) {
  set.seed(SEMILLA)   # gbm usa bag.fraction: fijar la semilla hace reproducible el ajuste
  m <- gbm::gbm(as.formula(paste("brote ~", paste(vars, collapse = " + "))),
                data = as.data.frame(tr), distribution = "bernoulli", n.trees = 250,
                interaction.depth = 3, shrinkage = 0.05, bag.fraction = 0.8,
                n.minobsinnode = 10, verbose = FALSE)
  predict(m, as.data.frame(te), n.trees = 250, type = "response")
}

# Bootstrap por bloques de semanas. Devuelve una matriz B x modelos con el AUC de
# cada réplica; en cada réplica se remuestrean semanas completas con reemplazo.
boot_auc <- function(te, scores, B) {
  filas <- split(seq_len(nrow(te)), te$week_start)
  y <- te$brote
  t(vapply(seq_len(B), function(b) {
    ix <- unlist(filas[sample(length(filas), replace = TRUE)], use.names = FALSE)
    vapply(scores, function(s) auc_rank(y[ix], s[ix]), numeric(1))
  }, numeric(length(scores))))
}

# Comparaciones de AUC: efecto del clima y ventaja sobre la persistencia
pares <- list(c("logistica_clima", "logistica_auto"),
              c("boosting_clima", "boosting_auto"),
              c("random_forest_clima", "random_forest_auto"),
              c("logistica_auto", "persistencia"),
              c("boosting_auto", "persistencia"))


# -----------------------------------------------------------------------------
# 3. Definición x horizonte
# -----------------------------------------------------------------------------
res_all <- list(); dif_all <- list()
for (def in DEFINICIONES) {
  depl <- etiquetar(dep, def)
  for (h in HORIZONTES) {
    ckpt <- file.path(DIR_SALIDA, sprintf("ckpt_S_%s_h%d.rds", def, h))
    if (REUSAR_CHECKPOINTS && file.exists(ckpt)) {
      cat(sprintf("[%s, h = %d] recuperado del punto de control\n", def, h))
      r <- readRDS(ckpt)
      res_all[[length(res_all) + 1]] <- r$res
      dif_all[[length(dif_all) + 1]] <- r$dif
      next
    }
    cat(sprintf("\nDefinición %s, horizonte %d\n", def, h))
    d <- prep_dep(depl, h)
    tr <- d |> filter(particion == "entrenamiento"); te <- d |> filter(particion == "prueba")
    te$departamento <- factor(te$departamento, levels = levels(tr$departamento))
    # Se alerta a la misma fracción de semanas que fueron brote en entrenamiento
    frac_alerta <- mean(tr$brote)
    cat(sprintf("Entrenamiento: %d (brotes %.1f%%) | Prueba: %d (brotes %.1f%%)\n",
                nrow(tr), 100 * mean(tr$brote), nrow(te), 100 * mean(te$brote)))

    scores <- list()
    for (con_clima in c(FALSE, TRUE)) {
      vars <- if (con_clima) vars_cli else vars_aut
      suf  <- if (con_clima) "_clima" else "_auto"
      m_lg <- suppressWarnings(glm(as.formula(paste("brote ~", paste(vars, collapse = " + "))),
                                   family = binomial(), data = tr))
      scores[[paste0("logistica", suf)]] <-
        suppressWarnings(predict(m_lg, newdata = te, type = "response"))
      m_rf <- ranger::ranger(as.formula(paste("brote_f ~", paste(vars, collapse = " + "))),
                             data = tr, probability = TRUE, num.trees = 400, min.node.size = 10,
                             seed = SEMILLA, num.threads = NTHREADS)
      scores[[paste0("random_forest", suf)]] <- predict(m_rf, data = te)$predictions[, "1"]
      scores[[paste0("boosting", suf)]] <- boost_brote(tr, te, vars)
    }
    # Referencia sin aprendizaje: persistencia (¿hubo brote hace h semanas?)
    scores[["persistencia"]] <- te$brote_l + 1e-6 * te$c_lag

    # Métricas puntuales
    y <- te$brote
    res <- bind_rows(lapply(names(scores), function(nm) {
      s <- scores[[nm]]
      alerta <- s >= quantile(s, 1 - frac_alerta)
      tibble(definicion = def, horizonte = h, modelo = nm, AUC = auc_rank(y, s),
             PR_AUC = pr_auc(y, s), Brier = mean((s - y)^2),
             precision_alerta = sum(alerta & y == 1) / sum(alerta),
             recall_alerta = sum(alerta & y == 1) / sum(y),
             lift = (sum(alerta & y == 1) / sum(alerta)) / mean(y),
             prev_entrenamiento = mean(tr$brote), prev_prueba = mean(y))
    }))

    # Intervalos de confianza por bootstrap de bloques de semanas
    cat("Bootstrap de bloques de semanas...\n")
    bt <- boot_auc(te, scores, B_BOOT)
    ic <- t(apply(bt, 2, quantile, probs = c(0.025, 0.975)))
    res <- res |> left_join(tibble(modelo = colnames(bt), AUC_lo = ic[, 1], AUC_hi = ic[, 2]),
                            by = "modelo")

    # Diferencias de AUC con su intervalo
    dif <- bind_rows(lapply(pares, function(p) {
      dd <- bt[, p[1]] - bt[, p[2]]
      tibble(definicion = def, horizonte = h, comparacion = paste(p[1], "-", p[2]),
             dAUC = auc_rank(y, scores[[p[1]]]) - auc_rank(y, scores[[p[2]]]),
             lo = unname(quantile(dd, 0.025)), hi = unname(quantile(dd, 0.975)),
             incluye_cero = lo <= 0 & hi >= 0)
    }))
    saveRDS(list(res = res, dif = dif), ckpt)
    res_all[[length(res_all) + 1]] <- res
    dif_all[[length(dif_all) + 1]] <- dif
  }
}

tabla <- bind_rows(res_all) |> mutate(across(where(is.numeric), ~ round(.x, 3)))
difs  <- bind_rows(dif_all) |> mutate(across(where(is.numeric), ~ round(.x, 3)))
write_csv(tabla, file.path(DIR_SALIDA, "tabla_sensibilidad_umbral.csv"))
write_csv(difs,  file.path(DIR_SALIDA, "tabla_diferencias_auc.csv"))

cat("\nAUC con IC 95 % (bootstrap de bloques de semanas)\n")
print(tabla |> select(definicion, horizonte, modelo, AUC, AUC_lo, AUC_hi, PR_AUC, lift, prev_prueba) |>
        arrange(definicion, horizonte, desc(AUC)), n = 200)
cat("\nDiferencias de AUC (incluye_cero = TRUE: sin evidencia de diferencia)\n")
print(difs, n = 100)


# -----------------------------------------------------------------------------
# 4. Figura
# -----------------------------------------------------------------------------
g10 <- tabla |>
  filter(modelo != "persistencia") |>
  mutate(clima = if_else(grepl("_clima$", modelo), "con clima", "sin clima"),
         familia = sub("_(auto|clima)$", "", modelo),
         definicion = factor(definicion, levels = DEFINICIONES,
                             labels = c("Percentil 75", "Percentil 90", "Media + 2 DE"))) |>
  ggplot(aes(factor(horizonte), AUC, color = familia, shape = clima)) +
  geom_hline(yintercept = 0.5, color = "grey60") +
  geom_pointrange(aes(ymin = AUC_lo, ymax = AUC_hi),
                  position = position_dodge(width = 0.7), size = 0.35) +
  facet_wrap(~ definicion) +
  labs(title = "Detección de brote: AUC según la definición de brote y la anticipación",
       subtitle = "Periodo de prueba 2022-2024. Barras: IC 95 % (bootstrap de bloques de semanas)",
       x = "Semanas de anticipación", y = "AUC", color = NULL, shape = NULL) + tema
mostrar(g10, "g10_sensibilidad_umbral.png", 11, 5)

cat("\nScript 04 finalizado. Tablas y figura en:", DIR_SALIDA, "\n")
