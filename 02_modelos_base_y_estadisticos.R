# =============================================================================
# Accidente ofídico en Colombia, 2015-2024
# Script 02. Líneas base y modelos estadísticos de conteo (horizonte de 1 semana)
#
# Requiere: 00_configuracion.R
# Entrada : dataset_modelado_brotes.rds (salida del script 01)
# Salida  : resultados_R/tabla_metricas_modelos.csv, coef_nb_clima_IRR.csv,
#           aporte_bloques_climaticos_AIC.csv, rezagos_precipitacion_IRR.csv,
#           rr_por_deciles_precipitacion.csv, predicciones_prueba.rds y
#           las figuras g05 y g06 (.png)
#
# Regla anti-fuga: para predecir la semana t solo se usa información hasta t-1.
# Las variables climáticas contemporáneas (z_*, precip_sum*w, etc.) se rezagan
# una semana en este script (sufijo _l1). Las columnas precip_mm, temp_* y
# rh_mean sin rezago no se utilizan.
#
# La precipitación entra en los modelos de cuatro formas:
#   (i)   anomalías acumuladas z_precip a 4, 13 y 26 semanas (efectos no lineales);
#   (ii)  lluvia semanal rezagada 1, 2, 4 y 8 semanas (rezagos distribuidos);
#   (iii) comparación explícita entre precipitación y temperatura más humedad
#         (sección 3c);
#   (iv)  razón entre casos observados y esperados por deciles de precipitación
#         (descriptivo).
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Configuración
# -----------------------------------------------------------------------------
paquetes <- c("dplyr", "readr", "tidyr", "ggplot2", "MASS", "mgcv", "scales")
faltan <- paquetes[!sapply(paquetes, requireNamespace, quietly = TRUE)]
if (length(faltan)) install.packages(faltan)
# MASS no se adjunta con library(): enmascararía dplyr::select. Se usa con MASS::.
for (p in setdiff(paquetes, "MASS")) library(p, character.only = TRUE)

source("00_configuracion.R")

# Dibuja los efectos suavizados de un GAM, sin el intercepto aleatorio por municipio
# ni las curvas estacionales propias de cada departamento (son 32 y no se pueden
# leer en un solo panel). El PNG se guarda primero, para que un fallo al dibujar en
# pantalla no impida conservar la figura.
graf_suaves <- function(m, archivo, ncol = 3) {
  etiquetas <- sapply(m$smooth, function(s) s$label)
  idx <- which(etiquetas != "s(municipio_id)" & !grepl("^s\\(semana_anio\\):", etiquetas))
  filas <- ceiling(length(idx) / ncol)
  dibuja <- function() {
    par(mfrow = c(filas, min(ncol, length(idx))))
    for (i in idx) plot(m, select = i, shade = TRUE, scale = 0,
                        ylab = "Efecto en log(tasa de incidencia)")
    par(mfrow = c(1, 1))
  }
  png(file.path(DIR_SALIDA, archivo), width = 2700, height = 900 * filas, res = 300)
  tryCatch(dibuja(), finally = dev.off())
  if (MOSTRAR_GRAFICOS) dibuja()
}

# Razón de tasas al duplicar la lluvia de la semana t-k, para los términos
# log1p(precip_mm_lagk) de un modelo. Como el término es log1p(lluvia),
# exp(coeficiente * log 2) es aproximadamente el cambio relativo de la tasa
# cuando la lluvia se duplica.
irr_rezagos <- function(m) {
  pt <- summary(m)$p.table
  rz <- pt[grepl("precip_mm_lag", rownames(pt)), , drop = FALSE]
  tibble(
    rezago_semanas = as.integer(sub(".*lag", "", sub("\\)", "", rownames(rz)))),
    razon_de_tasas = exp(rz[, "Estimate"] * log(2)),
    ic95_inf = exp((rz[, "Estimate"] - 1.96 * rz[, "Std. Error"]) * log(2)),
    ic95_sup = exp((rz[, "Estimate"] + 1.96 * rz[, "Std. Error"]) * log(2)),
    p_valor  = rz[, "Pr(>|z|)"]
  )
}

df <- readRDS(file.path(DIR_DATOS, "dataset_modelado_brotes.rds"))


# -----------------------------------------------------------------------------
# 1. Preparación de los datos
# -----------------------------------------------------------------------------
# 1.1 Rezago de una semana de las variables climáticas contemporáneas. El bloque 53
# (1 o 2 días) se excluyó del conjunto de datos, así que el rezago por posición
# equivale al bloque de 7 días anterior.
vars_contemp <- c("precip_sum4w", "precip_sum8w", "precip_sum13w", "precip_sum26w",
                  "rh_mean4w", "temp_mean4w",
                  "z_precip_4w", "z_precip_8w", "z_precip_13w", "z_precip_26w",
                  "z_rh_1w", "z_rh_4w", "z_temp_1w", "z_temp_4w")

df <- df |>
  arrange(municipio_id, week_start) |>
  group_by(municipio_id) |>
  mutate(across(all_of(vars_contemp), ~ dplyr::lag(.x, 1), .names = "{.col}_l1"),
         casos_lag1_chk = dplyr::lag(casos, 1)) |>
  ungroup() |>
  mutate(
    features_completas = as.logical(features_completas),
    sin2 = sin(4 * pi * semana_anio / 53),   # segundo armónico de la estacionalidad
    cos2 = cos(4 * pi * semana_anio / 53)
  )

# Control: casos_lag1 (calculado en Python) debe coincidir con el rezago por
# posición, salvo en la semana 1 de cada año, donde Python usó el bloque 53.
cat("Coincidencia entre casos_lag1 y el rezago por posición:",
    round(mean(df$casos_lag1 == df$casos_lag1_chk, na.rm = TRUE) * 100, 1), "%\n")

# 1.2 Municipios sin población válida. Se excluyen de los modelos con offset
# (Belén de Bajirá y Mapiripana: 45 casos de 51.343, es decir, 0,09 %), pero siguen
# contando en los agregados departamentales con que se define el brote.
excluidos <- df |>
  filter(is.na(log_poblacion_total) | poblacion_total <= 0) |>
  distinct(departamento, municipio)
cat("Municipios excluidos de los modelos municipales:\n")
print(excluidos)

preds_clima <- c("z_precip_4w_l1", "z_precip_13w_l1", "z_precip_26w_l1",
                 "z_temp_4w_l1", "z_rh_4w_l1")
preds_auto  <- c("casos_lag1", "casos_sum4w_prev", "casos_sum13w_prev",
                 "casos_dpto_sum4w_prev", "sin_semana", "cos_semana", "sin2", "cos2")
# Lluvia semanal rezagada, disponible hasta t-1: precip_mm_lagk es la lluvia de la semana t-k.
preds_lluvia <- c("precip_mm_lag1", "precip_mm_lag2", "precip_mm_lag4", "precip_mm_lag8")

dm <- df |>
  filter(features_completas, poblacion_total > 0, !is.na(log_poblacion_total)) |>
  drop_na(all_of(c(preds_auto, preds_clima, preds_lluvia, "casos_lag1_chk"))) |>
  mutate(departamento = factor(departamento), municipio_id = factor(municipio_id))

train <- dm |> filter(particion == "entrenamiento")
test  <- dm |> filter(particion == "prueba")
test$departamento <- factor(test$departamento, levels = levels(train$departamento))

# Un municipio con población válida solo en el periodo de prueba (Belén de Bajirá,
# COL.13.7_2, en 2024) no tiene historia de entrenamiento, así que sus efectos
# aleatorios no existen. Se retira de la evaluación para que todos los modelos se
# comparen sobre exactamente las mismas filas.
train <- droplevels(train)
n_antes <- nrow(test)
test  <- test |> filter(municipio_id %in% levels(train$municipio_id))
test$municipio_id <- factor(test$municipio_id, levels = levels(train$municipio_id))
cat("Filas de prueba retiradas por municipios sin entrenamiento:", n_antes - nrow(test), "\n")
cat(sprintf("\nEntrenamiento: %s filas | Prueba: %s filas\n",
            format(nrow(train), big.mark = "."), format(nrow(test), big.mark = ".")))


# -----------------------------------------------------------------------------
# 2. Líneas base (sin aprendizaje)
# -----------------------------------------------------------------------------
# B1, persistencia: lo ocurrido la semana anterior.
# B2, estacional por departamento: tasa histórica (casos por habitante) del
#     departamento en esa semana del año, multiplicada por la población.
# B3, media municipal: tasa histórica del municipio multiplicada por la población.
tasa_dep_sem <- train |>
  group_by(departamento, semana_anio) |>
  summarise(tasa_ds = sum(casos) / sum(poblacion_total), .groups = "drop")
tasa_mun <- train |>
  group_by(municipio_id) |>
  summarise(tasa_m = sum(casos) / sum(poblacion_total), .groups = "drop")

test <- test |>
  left_join(tasa_dep_sem, by = c("departamento", "semana_anio")) |>
  left_join(tasa_mun, by = "municipio_id") |>
  mutate(
    p_persistencia = casos_lag1,
    p_estac_depto  = coalesce(tasa_ds, 0) * poblacion_total,
    p_media_mun    = coalesce(tasa_m, 0)  * poblacion_total
  )


# -----------------------------------------------------------------------------
# 3. Modelos estadísticos de conteo (offset = log de la población)
# -----------------------------------------------------------------------------
f_auto <- casos ~ offset(log_poblacion_total) + departamento +
  sin_semana + cos_semana + sin2 + cos2 +
  log1p(casos_lag1) + log1p(casos_sum4w_prev) + log1p(casos_sum13w_prev) +
  log1p(casos_dpto_sum4w_prev)
f_clima <- update(f_auto, . ~ . + z_precip_4w_l1 + z_precip_13w_l1 + z_precip_26w_l1 +
                    z_temp_4w_l1 + z_rh_4w_l1)

cat("\nAjustando Poisson autorregresivo...\n")
m_pois <- glm(f_auto, family = poisson(), data = train)
disp <- sum(residuals(m_pois, type = "pearson")^2) / m_pois$df.residual
cat("Dispersión de Pearson del Poisson (1 = sin sobredispersión):", round(disp, 3), "\n")

cat("Ajustando binomial negativa autorregresiva...\n")
m_nb_auto <- MASS::glm.nb(f_auto, data = train)
cat("Ajustando binomial negativa autorregresiva con clima...\n")
m_nb_clima <- MASS::glm.nb(f_clima, data = train)
cat("theta de la binomial negativa (sin clima y con clima):",
    round(m_nb_auto$theta, 2), "y", round(m_nb_clima$theta, 2), "\n")

# ¿Aporta el clima? Prueba de razón de verosimilitud (modelos anidados) y AIC,
# ambos sobre el periodo de entrenamiento.
cat("\n¿El clima mejora el modelo? Prueba de razón de verosimilitud:\n")
print(anova(m_nb_auto, m_nb_clima))
cat("AIC sin clima:", round(AIC(m_nb_auto)), "| AIC con clima:", round(AIC(m_nb_clima)), "\n")

# Modelo aditivo con efectos no lineales del clima y estacionalidad cíclica
cat("\nAjustando el modelo aditivo binomial negativo (puede tardar unos minutos)...\n")
m_gam <- tryCatch(
  mgcv::bam(casos ~ offset(log_poblacion_total) + departamento +
              s(semana_anio, bs = "cc", k = 10) +
              log1p(casos_lag1) + log1p(casos_sum4w_prev) + log1p(casos_sum13w_prev) +
              log1p(casos_dpto_sum4w_prev) +
              s(z_precip_4w_l1, k = 6) + s(z_precip_13w_l1, k = 6) + s(z_precip_26w_l1, k = 6) +
              s(z_temp_4w_l1, k = 6) + s(z_rh_4w_l1, k = 6),
            family = mgcv::nb(), data = train, method = "fREML", discrete = TRUE),
  error = function(e) { message("El modelo aditivo falló: ", conditionMessage(e)); NULL })

# --- 3b. Modelos con intercepto aleatorio por municipio -----------------------
# La mayor parte de la heterogeneidad es entre municipios: uno de la Orinoquía o
# la Amazonía tiene una tasa base muy distinta a la de uno andino. El intercepto
# aleatorio la absorbe y evita que la binomial negativa sobreprediga. Se ajustan
# modelos con y sin clima para aislar el aporte climático.
cat("\nAjustando modelos con intercepto aleatorio por municipio (varios minutos)...\n")
base_re <- casos ~ offset(log_poblacion_total) +
  s(semana_anio, bs = "cc", k = 10) +
  log1p(casos_lag1) + log1p(casos_sum4w_prev) + log1p(casos_sum13w_prev) +
  log1p(casos_dpto_sum4w_prev) + s(municipio_id, bs = "re")

# Bloques de variables climáticas (todos rezagados: información hasta t-1)
b_precip_z <- ~ . + s(z_precip_4w_l1, k = 6) + s(z_precip_13w_l1, k = 6) + s(z_precip_26w_l1, k = 6)
b_precip_l <- ~ . + log1p(precip_mm_lag1) + log1p(precip_mm_lag2) +
  log1p(precip_mm_lag4) + log1p(precip_mm_lag8)
b_temp_rh  <- ~ . + s(z_temp_4w_l1, k = 6) + s(z_rh_4w_l1, k = 6)

clima_re  <- update(update(base_re, b_precip_z), b_temp_rh)   # todo el clima
precip_re <- update(update(base_re, b_precip_z), b_precip_l)  # solo precipitación (anomalías y rezagos)
temprh_re <- update(base_re, b_temp_rh)                       # solo temperatura y humedad

ajustar_bam <- function(f, fam) tryCatch(
  mgcv::bam(f, family = fam, data = train, method = "fREML", discrete = TRUE),
  error = function(e) { message("El ajuste falló: ", conditionMessage(e)); NULL })

m_re_auto     <- ajustar_bam(base_re,   poisson())
m_re_clima    <- ajustar_bam(clima_re,  poisson())
m_re_precip   <- ajustar_bam(precip_re, poisson())
m_re_temprh   <- ajustar_bam(temprh_re, poisson())
m_re_nb_clima <- ajustar_bam(clima_re,  mgcv::nb())

# --- 3c. ¿Qué variable climática aporta? Precipitación frente a temperatura y humedad
cat("\nAjustando modelos para separar anomalías y rezagos de lluvia...\n")
m_re_pz <- ajustar_bam(update(base_re, b_precip_z), poisson())   # solo anomalías acumuladas
m_re_pl <- ajustar_bam(update(base_re, b_precip_l), poisson())   # solo lluvia semanal rezagada

mods <- Filter(Negate(is.null), list(
  sin_clima = m_re_auto, temp_humedad = m_re_temprh,
  precip_solo_anomalias = m_re_pz, precip_solo_rezagos = m_re_pl,
  precip_anomalias_y_rezagos = m_re_precip, clima_completo_sin_rezagos = m_re_clima))
if ("sin_clima" %in% names(mods)) {
  aic_tab <- tibble(modelo = names(mods), AIC = sapply(mods, AIC),
                    edf = sapply(mods, function(m) sum(m$edf))) |>
    mutate(delta_AIC_vs_sin_clima = AIC - AIC[modelo == "sin_clima"])
  cat("\nAporte de cada bloque climático (menor AIC = mejor ajuste):\n")
  print(as.data.frame(aic_tab), digits = 8)
  write_csv(aic_tab, file.path(DIR_SALIDA, "aporte_bloques_climaticos_AIC.csv"))
}

# --- 3d. Robustez: ¿la señal de la lluvia sobrevive al control por año? -------
# Si la lluvia coincidiera con un cambio de nivel entre años (tendencia, cambios en
# la notificación, El Niño o La Niña), su efecto podría ser un efecto de periodo.
# Se repite la comparación con efectos fijos de año, solo con el entrenamiento.
cat("\nRobustez: modelos con efecto fijo de año...\n")
m_anio_base   <- ajustar_bam(update(base_re,   ~ . + factor(anio)), poisson())
m_anio_precip <- ajustar_bam(update(precip_re, ~ . + factor(anio)), poisson())
if (!is.null(m_anio_base) && !is.null(m_anio_precip)) {
  cat(sprintf("AIC controlando por año | sin lluvia: %.1f | con lluvia: %.1f | diferencia: %.1f\n",
              AIC(m_anio_base), AIC(m_anio_precip), AIC(m_anio_precip) - AIC(m_anio_base)))
  cat("Rezagos de lluvia (razón de tasas al duplicar), controlando por año:\n")
  print(irr_rezagos(m_anio_precip))
}

# --- 3e. Robustez: estacionalidad propia de cada departamento -----------------
# La suavización estacional del modelo base es común a todo el país. La lluvia
# semanal cruda sigue el ciclo estacional local (bimodal en los Andes, unimodal en
# la Orinoquía y la Amazonía), así que podría estar absorbiendo estacionalidad
# regional y no un efecto de la lluvia en sí. Aquí cada departamento tiene su
# propia curva estacional; si la lluvia sigue aportando, la señal es más sólida.
# Este bloque puede tardar varios minutos y se puede ejecutar por separado.
cat("\nRobustez: estacionalidad propia de cada departamento...\n")
base_ds <- casos ~ offset(log_poblacion_total) + departamento +
  s(semana_anio, by = departamento, bs = "cc", k = 8) +
  log1p(casos_lag1) + log1p(casos_sum4w_prev) + log1p(casos_sum13w_prev) +
  log1p(casos_dpto_sum4w_prev) + s(municipio_id, bs = "re")
m_ds_base   <- ajustar_bam(base_ds, poisson())
m_ds_precip <- ajustar_bam(update(base_ds, b_precip_l), poisson())
if (!is.null(m_ds_base) && !is.null(m_ds_precip)) {
  cat(sprintf("AIC con estacionalidad por departamento | sin lluvia: %.1f | con lluvia (rezagos): %.1f | diferencia: %.1f\n",
              AIC(m_ds_base), AIC(m_ds_precip), AIC(m_ds_precip) - AIC(m_ds_base)))
  cat("Rezagos de lluvia (razón de tasas al duplicar), con estacionalidad por departamento:\n")
  print(irr_rezagos(m_ds_precip))

  # Anomalías (desviaciones respecto de lo normal de cada lugar y época): con
  # estacionalidad propia por departamento, es la prueba más limpia del efecto de
  # la lluvia inusual, es decir, más o menos lluvia de la habitual.
  cat("\nAnomalías de lluvia, y de temperatura y humedad, con estacionalidad por departamento...\n")
  m_ds_pz    <- ajustar_bam(update(base_ds, b_precip_z), poisson())
  m_ds_clima <- ajustar_bam(update(update(base_ds, b_precip_z), b_temp_rh), poisson())
  if (!is.null(m_ds_pz))
    cat(sprintf("Diferencia de AIC con anomalías de lluvia: %.1f\n", AIC(m_ds_pz) - AIC(m_ds_base)))
  if (!is.null(m_ds_clima)) {
    cat(sprintf("Diferencia de AIC con anomalías de lluvia, temperatura y humedad: %.1f\n",
                AIC(m_ds_clima) - AIC(m_ds_base)))
    st3 <- summary(m_ds_clima)$s.table
    print(st3[grepl("z_", rownames(st3)), , drop = FALSE])
    tryCatch(graf_suaves(m_ds_clima, "g06d_efectos_anomalias_estacionalidad_depto.png"),
             error = function(e) message("No se pudo graficar: ", conditionMessage(e)))
  }
}

if (!is.null(m_re_precip)) {
  cat("\nEfectos suavizados de la precipitación (grados de libertad y valor p):\n")
  st <- summary(m_re_precip)$s.table
  print(st[grepl("precip", rownames(st)), , drop = FALSE])

  # Curvas de respuesta a la precipitación y a la estacionalidad, con banda de confianza
  tryCatch(graf_suaves(m_re_precip, "g06_efectos_precipitacion.png"),
           error = function(e) message("No se pudo graficar: ", conditionMessage(e)))

  # Rezagos distribuidos: cambio de la tasa al duplicar la lluvia de la semana t-k
  lag_tab <- irr_rezagos(m_re_precip)
  cat("\nRezagos de la lluvia semanal (razón de tasas al duplicar la lluvia):\n")
  print(lag_tab)
  write_csv(lag_tab, file.path(DIR_SALIDA, "rezagos_precipitacion_IRR.csv"))
  g6c <- ggplot(lag_tab, aes(factor(rezago_semanas), razon_de_tasas)) +
    geom_hline(yintercept = 1, color = "grey50") +
    geom_pointrange(aes(ymin = ic95_inf, ymax = ic95_sup), color = "#2b6a99") +
    labs(title = "Efecto de la lluvia semanal rezagada sobre la incidencia",
         subtitle = "Razón de tasas al duplicar (aprox.) la lluvia de la semana t-k; ajustado por estacionalidad, municipio y casos previos",
         x = "Rezago k (semanas antes)", y = "Razón de tasas (IC 95 %)") +
    theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"),
                                          plot.subtitle = element_text(size = 9))
  mostrar(g6c, "g06c_rezagos_precipitacion.png", 8, 4.5)
}

if (!is.null(m_re_clima)) {
  cat("\nEfectos suavizados del clima completo (grados de libertad y valor p):\n")
  print(summary(m_re_clima)$s.table)
  tryCatch(graf_suaves(m_re_clima, "g06_efectos_clima_completo.png"),
           error = function(e) message("No se pudo graficar: ", conditionMessage(e)))
}

# --- Descriptivo: razón entre casos observados y esperados por deciles de precipitación
# El esperado es la tasa media del propio municipio por su población, lo que controla
# el nivel base de cada municipio. RR > 1 indica más casos de los habituales para ese
# municipio. Es un análisis descriptivo que usa todos los años (entrenamiento y prueba).
rr_precip <- bind_rows(lapply(
  c(`Anomalía 4 sem.` = "z_precip_4w_l1", `Anomalía 13 sem.` = "z_precip_13w_l1",
    `Anomalía 26 sem.` = "z_precip_26w_l1"),
  function(v) {
    dm |> group_by(municipio_id) |>
      mutate(esperado = sum(casos) / sum(poblacion_total) * poblacion_total) |> ungroup() |>
      mutate(decil = ntile(.data[[v]], 10)) |> group_by(decil) |>
      summarise(z_medio = mean(.data[[v]]), casos = sum(casos), esperado = sum(esperado),
                .groups = "drop") |>
      mutate(RR = casos / esperado, lo = exp(log(RR) - 1.96 / sqrt(casos)),
             hi = exp(log(RR) + 1.96 / sqrt(casos)))
  }), .id = "ventana") |>
  mutate(ventana = factor(ventana, levels = c("Anomalía 4 sem.", "Anomalía 13 sem.", "Anomalía 26 sem.")))
write_csv(rr_precip, file.path(DIR_SALIDA, "rr_por_deciles_precipitacion.csv"))
g6b <- ggplot(rr_precip, aes(z_medio, RR)) +
  geom_hline(yintercept = 1, color = "grey50") +
  geom_pointrange(aes(ymin = lo, ymax = hi), color = "#2b6a99") +
  geom_line(color = "#2b6a99", alpha = 0.5) +
  facet_wrap(~ ventana, scales = "free_x") +
  labs(title = "Casos observados y esperados según la anomalía de precipitación previa",
       subtitle = "Cada punto es un decil (todos los años). Esperado = tasa media del municipio x población.",
       x = "Anomalía de precipitación (desviaciones estándar respecto de la climatología)",
       y = "Razón observado / esperado (IC 95 %)") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"),
                                        plot.subtitle = element_text(size = 9))
mostrar(g6b, "g06b_razon_obs_esperado_precipitacion.png", 11, 4.5)

# Predicciones en el periodo de prueba, en escala de conteo esperado
test$p_pois     <- predict(m_pois,     newdata = test, type = "response")
test$p_nb_auto  <- predict(m_nb_auto,  newdata = test, type = "response")
test$p_nb_clima <- predict(m_nb_clima, newdata = test, type = "response")
if (!is.null(m_gam))         test$p_gam_clima   <- as.numeric(predict(m_gam,         newdata = test, type = "response"))
if (!is.null(m_re_auto))     test$p_re_auto     <- as.numeric(predict(m_re_auto,     newdata = test, type = "response"))
if (!is.null(m_re_clima))    test$p_re_clima    <- as.numeric(predict(m_re_clima,    newdata = test, type = "response"))
if (!is.null(m_re_nb_clima)) test$p_re_nb_clima <- as.numeric(predict(m_re_nb_clima, newdata = test, type = "response"))
if (!is.null(m_re_precip))   test$p_re_precip   <- as.numeric(predict(m_re_precip,   newdata = test, type = "response"))
if (!is.null(m_re_temprh))   test$p_re_temprh   <- as.numeric(predict(m_re_temprh,   newdata = test, type = "response"))

# Coeficientes de la binomial negativa con clima como razones de tasas (IRR),
# con intervalos de confianza de Wald al 95 %.
co <- summary(m_nb_clima)$coefficients
irr <- tibble(
  termino  = rownames(co),
  IRR      = exp(co[, "Estimate"]),
  IC95_inf = exp(co[, "Estimate"] - 1.96 * co[, "Std. Error"]),
  IC95_sup = exp(co[, "Estimate"] + 1.96 * co[, "Std. Error"]),
  p_valor  = co[, "Pr(>|z|)"]
) |> filter(!grepl("^departamento", termino))
cat("\nIRR de la binomial negativa con clima (sin los efectos de departamento):\n")
print(irr, n = 30)
write_csv(irr, file.path(DIR_SALIDA, "coef_nb_clima_IRR.csv"))


# -----------------------------------------------------------------------------
# 4. Evaluación en el periodo de prueba (2022-2024)
# -----------------------------------------------------------------------------
# Agrega casos observados y predichos por las claves indicadas (municipio,
# departamento o nación), conservando el umbral y la etiqueta de brote.
agrega <- function(d, col, claves) {
  d |> group_by(across(all_of(claves))) |>
    summarise(obs = sum(casos), pred = sum(.data[[col]]),
              umbral = first(umbral_brote_dpto), brote = first(brote_dpto),
              .groups = "drop")
}

resumen <- function(a, nivel, modelo) {
  tibble(modelo = modelo, nivel = nivel,
         MAE  = mean(abs(a$obs - a$pred)),
         RMSE = sqrt(mean((a$obs - a$pred)^2)),
         R2   = 1 - sum((a$obs - a$pred)^2) / sum((a$obs - mean(a$obs))^2),
         dev_poisson = pois_dev(a$obs, a$pred),
         casos_obs = sum(a$obs), casos_pred = round(sum(a$pred)))
}

modelos <- c(persistencia = "p_persistencia", estacional_depto = "p_estac_depto",
             media_municipal = "p_media_mun", poisson_auto = "p_pois",
             nb_auto = "p_nb_auto", nb_clima = "p_nb_clima",
             gam_clima = "p_gam_clima", pois_re_auto = "p_re_auto",
             pois_re_clima = "p_re_clima", nb_re_clima = "p_re_nb_clima",
             pois_re_precip = "p_re_precip", pois_re_temprh = "p_re_temprh")
modelos <- modelos[modelos %in% names(test)]

tabla <- bind_rows(lapply(names(modelos), function(nm) {
  col <- modelos[[nm]]
  bind_rows(
    resumen(agrega(test, col, c("municipio_id", "week_start")), "municipio-semana", nm),
    resumen(agrega(test, col, c("departamento", "week_start")), "departamento-semana", nm),
    resumen(agrega(test, col, "week_start"), "nacional-semana", nm)
  )
}))

# AUC de la detección de brote, usando como puntaje la predicción departamental
# dividida por el umbral del canal endémico.
auc_tab <- bind_rows(lapply(names(modelos), function(nm) {
  a <- agrega(test, modelos[[nm]], c("departamento", "week_start")) |> drop_na(umbral, brote)
  tibble(modelo = nm, AUC_brote_depto = auc_rank(a$brote, a$pred / a$umbral))
}))

tabla_final <- tabla |> left_join(auc_tab, by = "modelo") |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))
cat("\nMétricas en el periodo de prueba:\n")
print(tabla_final |> arrange(nivel, dev_poisson), n = 50)
write_csv(tabla_final, file.path(DIR_SALIDA, "tabla_metricas_modelos.csv"))


# -----------------------------------------------------------------------------
# 5. Gráfico: casos observados y predichos (nacional, semana a semana, 2022-2024)
# -----------------------------------------------------------------------------
sel <- intersect(c("persistencia", "poisson_auto", "pois_re_auto", "pois_re_precip"), names(modelos))
nac <- bind_rows(
  agrega(test, "p_persistencia", "week_start") |> transmute(week_start, valor = obs, serie = "Observado"),
  bind_rows(lapply(sel, function(nm)
    agrega(test, modelos[[nm]], "week_start") |> transmute(week_start, valor = pred, serie = nm)))
)
g5 <- ggplot(nac, aes(week_start, valor, color = serie, linewidth = serie == "Observado")) +
  geom_line() + scale_linewidth_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.4), guide = "none") +
  scale_color_manual(values = c(Observado = "black", persistencia = "grey60",
                                poisson_auto = "#2b6a99", pois_re_auto = "#27ae60", pois_re_precip = "#c0392b"),
                     name = NULL) +
  labs(title = "Casos nacionales por semana en el periodo de prueba",
       subtitle = "Predicción a 1 semana con información hasta t-1", x = NULL, y = "Casos") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
mostrar(g5, "g05_obs_vs_pred.png", 9, 4.5)

saveRDS(test |> select(municipio_id, departamento, week_start, casos, brote_dpto,
                       umbral_brote_dpto, starts_with("p_")),
        file.path(DIR_SALIDA, "predicciones_prueba.rds"))
cat("\nScript 02 finalizado. Siguiente paso: script 03 (aprendizaje automático y clasificación de brotes).\n")
