# =============================================================================
# Accidente ofídico en Colombia, 2015-2024
# Script 00. Configuración compartida
#
# Los scripts 01 a 04 cargan este archivo con source("00_configuracion.R").
# Reúne las rutas, los parámetros del análisis, las semillas y las funciones
# que se repetían en varios scripts. Para cambiar un parámetro (por ejemplo, el
# mínimo de casos para declarar un brote) basta con editarlo aquí.
#
# Los scripts se ejecutan con el directorio de trabajo en la carpeta del
# repositorio (en RStudio: Session > Set Working Directory > To Source File
# Location, o abrir el proyecto .Rproj).
# =============================================================================

options(width = 250, pillar.width = Inf, tibble.width = Inf, scipen = 999)

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
suppressPackageStartupMessages(library(ggplot2))


# -----------------------------------------------------------------------------
# Rutas
# -----------------------------------------------------------------------------
# Carpeta con dataset_modelado.csv y los archivos intermedios (por defecto, el
# directorio de trabajo). Los resultados se guardan en la subcarpeta resultados_R.
DIR_DATOS  <- "."
DIR_SALIDA <- file.path(DIR_DATOS, "resultados_R")
dir.create(DIR_SALIDA, showWarnings = FALSE, recursive = TRUE)


# -----------------------------------------------------------------------------
# Parámetros del análisis
# -----------------------------------------------------------------------------
ANIO_FIN_ENTRENAMIENTO <- 2021                          # coincide con la columna `particion`
ANIOS_BASE_CANAL       <- 2015:ANIO_FIN_ENTRENAMIENTO   # años con que se construye el canal endémico
ANIOS_EXCLUIDOS_CANAL  <- 2020                          # pandemia de COVID-19: subnotificación atípica
VENTANA_SEMANAS        <- 2                             # el umbral de la semana s usa las semanas s-2 a s+2
PERCENTIL_BROTE        <- 0.75                          # límite superior de la zona esperada del canal
MIN_CASOS_BROTE        <- 3                             # evita declarar brote con 1 o 2 casos en departamentos pequeños
HORIZONTES             <- c(1, 4, 8)                    # semanas de anticipación


# -----------------------------------------------------------------------------
# Ejecución y reproducibilidad
# -----------------------------------------------------------------------------
SEMILLA            <- 2024      # modelos (gbm, ranger) y bootstrap
MOSTRAR_GRAFICOS   <- TRUE      # TRUE dibuja cada gráfico en RStudio además de guardarlo
REUSAR_CHECKPOINTS <- TRUE      # TRUE retoma los puntos de control ckpt_*.rds si existen
NTHREADS           <- max(1, min(4, parallel::detectCores() - 1))
set.seed(SEMILLA)


# -----------------------------------------------------------------------------
# Funciones comunes
# -----------------------------------------------------------------------------
tema <- theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))

# Guarda un gráfico de ggplot en PNG (300 dpi) y, si MOSTRAR_GRAFICOS es TRUE, lo
# dibuja en RStudio. print() es necesario porque, con source() o dentro de una
# función, ggplot no dibuja el gráfico por sí solo.
mostrar <- function(g, archivo, w = 9, h = 4.5) {
  ggsave(file.path(DIR_SALIDA, archivo), g, width = w, height = h, dpi = 300)
  if (MOSTRAR_GRAFICOS) print(g)
  invisible(g)
}

# Desviación de Poisson media (menor es mejor)
pois_dev <- function(y, mu) {
  mu <- pmax(mu, 1e-8)
  2 * mean(ifelse(y > 0, y * log(y / mu), 0) - (y - mu))
}

# AUC por rangos (equivalente a la U de Mann-Whitney), sin paquetes adicionales
auc_rank <- function(y, s) {
  r <- rank(s); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Instala, si falta, un paquete en su versión binaria en Windows (evita compilar)
# y devuelve TRUE si queda disponible.
instalar_bin <- function(p) {
  if (!requireNamespace(p, quietly = TRUE))
    try(suppressWarnings(install.packages(
      p, quiet = TRUE,
      type = if (.Platform$OS.type == "windows") "binary" else getOption("pkgType"))), silent = TRUE)
  requireNamespace(p, quietly = TRUE)
}
