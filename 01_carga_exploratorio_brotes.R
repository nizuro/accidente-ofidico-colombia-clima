# =============================================================================
# Accidente ofídico en Colombia, 2015-2024
# Script 01. Carga, verificación, análisis exploratorio y definición de brote
#
# Requiere: 00_configuracion.R
# Entrada : dataset_modelado.csv (generado en el cuaderno de Colab, sección 8)
# Salida  : resultados_R/g01 a g04 (.png), canal_endemico_departamentos.csv,
#           departamento_semana_brotes.csv y dataset_modelado_brotes.rds
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Configuración (parámetros y rutas en 00_configuracion.R)
# -----------------------------------------------------------------------------
paquetes <- c("dplyr", "readr", "tidyr", "ggplot2", "lubridate", "stringr", "scales")
faltan <- paquetes[!sapply(paquetes, requireNamespace, quietly = TRUE)]
if (length(faltan)) install.packages(faltan)
invisible(lapply(paquetes, library, character.only = TRUE))

source("00_configuracion.R")
tema <- tema + theme(panel.grid.minor = element_blank())


# -----------------------------------------------------------------------------
# 1. Carga y verificación
# -----------------------------------------------------------------------------
df <- read_csv(
  file.path(DIR_DATOS, "dataset_modelado.csv"),
  col_types = cols(
    municipio_id = col_character(),
    dane_code    = col_character(),   # como texto, para no perder los ceros iniciales
    week_start   = col_date(format = "%Y-%m-%d"),
    .default     = col_guess()
  )
) |>
  mutate(
    dane_code = str_pad(dane_code, 5, pad = "0"),
    particion = factor(particion, levels = c("entrenamiento", "prueba"))
  )

cat("Dimensiones:", dim(df), "\n")
stopifnot(!anyDuplicated(df[, c("municipio_id", "week_start")]))
stopifnot(all(df$casos >= 0), all(df$casos == round(df$casos)))

# Partición temporal: debe ser 2015-2021 (entrenamiento) y 2022-2024 (prueba).
df |>
  group_by(particion) |>
  summarise(desde = min(anio), hasta = max(anio), filas = n(),
            municipios = n_distinct(municipio_id), casos = sum(casos),
            .groups = "drop") |>
  print()

# Valores faltantes, solo de las variables que tienen alguno.
faltantes <- df |>
  summarise(across(everything(), ~ sum(is.na(.x)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_na") |>
  filter(n_na > 0)
print(faltantes)

# Continuidad semanal. Las semanas son bloques de 7 días desde el 1 de enero y
# se excluyó el bloque 53 (dura 1 o 2 días), de modo que al cambiar de año el
# salto entre filas consecutivas es de 8 o 9 días. Este detalle importa en los
# scripts siguientes, porque un rezago de una fila equivale entonces a 8 o 9 días.
cat("\nFilas por semana del año (la 53 debería tener pocas o ninguna):\n")
print(table(df$semana_anio))

cat("\nTamaño de los saltos entre filas consecutivas de un mismo municipio (días):\n")
df |>
  arrange(municipio_id, week_start) |>
  group_by(municipio_id) |>
  mutate(salto_dias = as.integer(week_start - lag(week_start))) |>
  ungroup() |>
  filter(!is.na(salto_dias)) |>
  count(salto_dias) |>
  print()

# Municipios sin población válida (ausente o no positiva). Sus filas no entran
# en los modelos municipales, que usan el logaritmo de la población como offset.
cat("\nPoblación faltante o no positiva (municipio y años afectados):\n")
df |>
  filter(is.na(poblacion_total) | poblacion_total <= 0 | is.na(log_poblacion_total)) |>
  group_by(departamento, municipio) |>
  summarise(filas = n(), anios = paste(sort(unique(anio)), collapse = ","),
            pobl_min = suppressWarnings(min(poblacion_total, na.rm = TRUE)),
            casos = sum(casos), .groups = "drop") |>
  arrange(desc(filas)) |>
  print(n = 30)


# -----------------------------------------------------------------------------
# 2. Análisis exploratorio
# -----------------------------------------------------------------------------
# 2.1 Serie nacional semanal
nacional <- df |>
  group_by(week_start) |>
  summarise(casos = sum(casos), .groups = "drop")

g1 <- ggplot(nacional, aes(week_start, casos)) +
  geom_line(color = "#2b6a99", linewidth = 0.4) +
  geom_smooth(span = 0.15, se = FALSE, color = "#c0392b", linewidth = 0.7) +
  labs(title = "Accidente ofídico en Colombia: casos semanales",
       subtitle = "Línea roja: tendencia suavizada",
       x = NULL, y = "Casos por semana") + tema
mostrar(g1, "g01_serie_nacional.png", 9, 4.5)

# 2.2 Estacionalidad: casos nacionales por semana del año, una línea por año
estacionalidad <- df |>
  group_by(anio, semana_anio) |>
  summarise(casos = sum(casos), .groups = "drop")

g2 <- ggplot(estacionalidad, aes(semana_anio, casos, group = anio, color = factor(anio))) +
  geom_line(alpha = 0.8) +
  scale_color_viridis_d(name = "Año") +
  labs(title = "Estacionalidad nacional por semana del año",
       x = "Semana del año", y = "Casos") + tema
mostrar(g2, "g02_estacionalidad.png", 9, 4.5)

# 2.3 Los 15 departamentos con más casos
dep_tot <- df |>
  group_by(departamento) |>
  summarise(casos = sum(casos), .groups = "drop") |>
  slice_max(casos, n = 15)

g3 <- ggplot(dep_tot, aes(reorder(departamento, casos), casos)) +
  geom_col(fill = "#2b6a99") + coord_flip() +
  scale_y_continuous(labels = label_number(big.mark = ".")) +
  labs(title = "Departamentos con más casos (2015-2024)", x = NULL, y = "Casos") + tema
mostrar(g3, "g03_top_departamentos.png", 7, 5)

# 2.4 Distribución de la variable respuesta en el nivel municipio-semana
cat("\nProporción de municipio-semanas con 0 casos:",
    round(mean(df$casos == 0) * 100, 1), "%\n")
cat("Media y varianza de los casos:", round(mean(df$casos), 3), "y", round(var(df$casos), 3), "\n")
# La razón varianza/media marginal es solo un indicio de sobredispersión, porque
# la media varía mucho entre municipios. La sobredispersión propiamente dicha se
# evalúa en el script 02, con los residuos de un modelo con población y covariables.

# 2.5 Relación entre el clima y los casos con rezago, a nivel nacional semanal
clima_nac <- df |>
  group_by(week_start) |>
  summarise(casos = sum(casos),
            z_precip_13w = mean(z_precip_13w, na.rm = TRUE),
            z_temp_4w    = mean(z_temp_4w,    na.rm = TRUE),
            z_rh_4w      = mean(z_rh_4w,      na.rm = TRUE), .groups = "drop") |>
  arrange(week_start)

correlacion_cruzada <- sapply(c("z_precip_13w", "z_temp_4w", "z_rh_4w"), function(v) {
  r <- ccf(clima_nac[[v]], clima_nac$casos, lag.max = 12, plot = FALSE, na.action = na.omit)
  idx <- r$lag[, 1, 1] <= 0   # rezago negativo: el clima antecede a los casos
  setNames(r$acf[idx, 1, 1], r$lag[idx, 1, 1])
})
cat("\nCorrelación cruzada (el clima antecede a los casos; columnas = rezago en semanas):\n")
print(round(t(correlacion_cruzada), 2))


# -----------------------------------------------------------------------------
# 3. Definición de brote con un canal endémico (Bortman, 1999)
# -----------------------------------------------------------------------------
# Se trabaja en el nivel departamento-semana: a escala municipal casi todas las
# semanas tienen cero casos y no es posible estimar un umbral estable. El umbral
# se calcula solo con los años base (2015-2021 sin 2020), por lo que no hay fuga
# de información hacia el periodo de prueba.

dep_sem <- df |>
  group_by(departamento, anio, semana_anio, week_start) |>
  summarise(casos = sum(casos), .groups = "drop")

base <- dep_sem |>
  filter(anio %in% ANIOS_BASE_CANAL, !anio %in% ANIOS_EXCLUIDOS_CANAL)

n_anios_base <- n_distinct(base$anio)
cat("\nAños en la base del canal:", n_anios_base,
    "(Bortman recomienda entre 5 y 11):", paste(sort(unique(base$anio)), collapse = ", "), "\n")
stopifnot(n_anios_base >= 5)

# Para cada semana s se toman los casos de las semanas s-2 a s+2 de todos los años
# base y se calcula su mediana, su percentil 75 y su media. La ventana es circular
# (módulo 53) para evitar bordes irregulares en las primeras y últimas semanas.
semanas <- 1:53
umbral_semana <- function(dat, s) {
  ventana <- ((s - 1 + (-VENTANA_SEMANAS:VENTANA_SEMANAS)) %% 53) + 1
  x <- dat$casos[dat$semana_anio %in% ventana]
  c(mediana = median(x),
    q75     = unname(quantile(x, PERCENTIL_BROTE, type = 7)),
    media   = mean(x))
}

canal <- base |>
  group_by(departamento) |>
  group_modify(~ {
    resumen <- t(sapply(semanas, function(s) umbral_semana(.x, s)))
    tibble(semana_anio = semanas, as_tibble(resumen))
  }) |>
  ungroup()

# Etiqueta de brote para todo el periodo, entrenamiento y prueba. El umbral es el
# percentil 75, pero nunca menor que MIN_CASOS_BROTE - 1: se declara brote cuando
# los casos superan el umbral, es decir, cuando son al menos MIN_CASOS_BROTE.
dep_sem <- dep_sem |>
  left_join(canal, by = c("departamento", "semana_anio")) |>
  mutate(
    umbral_brote = pmax(q75, MIN_CASOS_BROTE - 1),
    brote        = as.integer(casos > umbral_brote)
  )

cat("\nSemanas-departamento en brote, por periodo:\n")
dep_sem |>
  mutate(periodo = if_else(anio <= ANIO_FIN_ENTRENAMIENTO, "entrenamiento", "prueba")) |>
  group_by(periodo) |>
  summarise(n = n(), brotes = sum(brote), prevalencia = round(mean(brote) * 100, 1),
            .groups = "drop") |>
  print()

# Canal endémico de los cuatro departamentos con más casos, en tres años de ejemplo
deptos_demo <- dep_tot$departamento[1:4]
g4 <- dep_sem |>
  filter(departamento %in% deptos_demo, anio %in% c(2019, 2022, 2023)) |>
  ggplot(aes(semana_anio)) +
  geom_ribbon(aes(ymin = 0, ymax = umbral_brote), fill = "grey85") +
  geom_line(aes(y = mediana), color = "grey40", linetype = "dashed") +
  geom_line(aes(y = casos), color = "#c0392b") +
  geom_point(data = ~ filter(.x, brote == 1), aes(y = casos), color = "#c0392b", size = 1.4) +
  facet_grid(departamento ~ anio, scales = "free_y") +
  labs(title = "Canal endémico por departamento",
       subtitle = "Zona gris: casos esperados (hasta el percentil 75). Puntos: semanas en brote",
       x = "Semana del año", y = "Casos") + tema
mostrar(g4, "g04_canal_endemico.png", 10, 7)


# -----------------------------------------------------------------------------
# 4. Guardado
# -----------------------------------------------------------------------------
write_csv(canal,   file.path(DIR_SALIDA, "canal_endemico_departamentos.csv"))
write_csv(dep_sem, file.path(DIR_SALIDA, "departamento_semana_brotes.csv"))

df <- df |>
  left_join(dep_sem |> select(departamento, week_start, casos_dpto = casos,
                              umbral_brote_dpto = umbral_brote, brote_dpto = brote),
            by = c("departamento", "week_start"))
saveRDS(df, file.path(DIR_DATOS, "dataset_modelado_brotes.rds"))

cat("\nScript 01 finalizado. Siguiente paso: script 02 (modelos base y estadísticos).\n")
