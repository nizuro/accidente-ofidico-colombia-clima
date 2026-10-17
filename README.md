# Accidente ofídico en Colombia, 2015-2024: clima, brotes y predicción

Código del Trabajo Fin de Máster del Máster en Big Data y Ciencia de Datos (Universidad Internacional de Valencia, VIU). El trabajo analiza la relación entre la precipitación y el accidente ofídico (evento 100 del Sivigila) en los municipios de Colombia. También evalúa si el clima mejora la predicción del conteo semanal de casos y de la detección de brotes departamentales a 1, 4 y 8 semanas.

**Autor:** Nicolás Zúñiga · **Año:** 2026

## Contenido del repositorio

| Archivo | Descripción |
|---|---|
| `01_Notebook_Colab_construccion_datos.ipynb` | Construcción del panel municipio-semana (casos, clima, población), análisis espacial (suavizado bayesiano empírico, I de Moran, LISA, Gi*) y exportación de `dataset_modelado.csv`. Se ejecuta en Google Colab (Python). |
| `00_configuracion.R` | Rutas, parámetros del análisis, semillas y funciones comunes. Lo cargan los scripts 01 a 04. |
| `01_carga_exploratorio_brotes.R` | Carga y verificación del panel, análisis exploratorio y definición de brote mediante canal endémico por departamento. |
| `02_modelos_base_y_estadisticos.R` | Líneas base y modelos estadísticos de conteo (Poisson, binomial negativa, modelos aditivos) con y sin clima, horizonte de 1 semana. |
| `03_aprendizaje_automatico_horizontes.R` | Boosting de Poisson para el conteo; regresión logística, bosque aleatorio y boosting para el brote, a 1, 4 y 8 semanas. |
| `04_sensibilidad_umbral_e_intervalos.R` | Sensibilidad a la definición de brote (percentil 75, percentil 90, media + 2 desviaciones estándar) e intervalos de confianza por bootstrap de bloques de semanas. |
| `transparencia_ia/prompts_utilizados.md` | Registro de los prompts utilizados con la herramienta de inteligencia artificial generativa. |

## Uso de inteligencia artificial generativa

El trabajo contó con el apoyo de Claude (Anthropic, 2026) en la redacción preliminar, la escritura y revisión del código y la comprobación de referencias. El registro de los prompts utilizados, con fecha y propósito, está en [`transparencia_ia/prompts_utilizados.md`](transparencia_ia/prompts_utilizados.md). Todos los resultados provienen de la ejecución del código sobre los datos originales.

## Datos

Los datos crudos no se incluyen en el repositorio por su tamaño y licencia. Deben descargarse de las fuentes originales y colocarse como indica el cuaderno de Colab.

| Dato | Fuente y enlace | Condiciones de uso |
|---|---|---|
| Casos de accidente ofídico (evento 100), 2015-2024 | Instituto Nacional de Salud (INS), Sivigila. Portal de [Datos Abiertos del INS](https://www.ins.gov.co/Transparencia/informacion-de-interes/Paginas/Datos-Abiertos.aspx), con los conjuntos publicados en [datos.gov.co](https://www.datos.gov.co/). Se usan los archivos anuales `Datos_AAAA_Accidente_ofidico.xlsx`. | Datos abiertos (Ley 1712 de 2014). Revisar la licencia de cada conjunto en datos.gov.co. |
| Precipitación diaria, 0,05° | [CHIRPS v2.0](https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/tifs/p05/) (Climate Hazards Center, UC Santa Barbara). Referencia: Funk et al. (2015), <https://doi.org/10.1038/sdata.2015.66>. | Dominio público. |
| Temperatura y humedad (temperatura del aire y del punto de rocío) | [ERA5-Land hourly data from 1950 to present](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land), Copernicus Climate Change Service, <https://doi.org/10.24381/cds.e2161bac>. Referencia: Muñoz-Sabater et al. (2021), <https://doi.org/10.5194/essd-13-4349-2021>. Requiere una cuenta y una clave de la API del Climate Data Store. | CC BY 4.0. |
| Población municipal total y rural | DANE, [Proyecciones de población](https://www.dane.gov.co/index.php/estadisticas-por-tema/demografia-y-poblacion/proyecciones-de-poblacion): proyección municipal por área 2018-2042 ([PPED-AreaMun-2018-2042_VP.xlsx](https://www.dane.gov.co/files/censo2018/proyecciones-de-poblacion/Municipal/PPED-AreaMun-2018-2042_VP.xlsx)) y retroproyección 2005-2017 (`DCD-area-proypoblacion-Mun-2005-2017_VP.xlsx`, en la misma página). | [Términos y condiciones del DANE](https://www.dane.gov.co/index.php/terrminos-y-condiciones-dane). |
| Límites municipales | [GADM 4.1, nivel 2](https://gadm.org/download_country.html) (Colombia). El cuaderno lo descarga de <https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_COL_2.json>. | Libre para uso académico y no comercial; la redistribución requiere permiso. Por eso el repositorio no incluye los polígonos. |

Todos los enlaces se consultaron el 17 de octubre de 2026. CHIRPS v2.0 dejará de producirse después de diciembre de 2026; los archivos ya publicados siguen disponibles.

## Cómo reproducir los resultados

1. Ejecutar el cuaderno de Colab de principio a fin. Genera `dataset_modelado.csv` (panel de 578 088 filas: municipio por semana, 2015-2024).
2. Descargar `dataset_modelado.csv` y copiarlo en la carpeta del repositorio.
3. Abrir RStudio, fijar el directorio de trabajo en la carpeta del repositorio y ejecutar en orden:

```r
source("01_carga_exploratorio_brotes.R")
source("02_modelos_base_y_estadisticos.R")
source("03_aprendizaje_automatico_horizontes.R")
source("04_sensibilidad_umbral_e_intervalos.R")
```

Los resultados (tablas CSV y figuras PNG) se guardan en `resultados_R/`. Cada script instala los paquetes que le faltan.

**Tiempos y puntos de control.** Los scripts 03 y 04 son los más lentos (el 04 con 1 000 réplicas bootstrap). Guardan un punto de control por horizonte (`ckpt_*.rds`) y, si se interrumpen, retoman donde quedaron. Para recalcular desde cero, borrar esos archivos o poner `REUSAR_CHECKPOINTS <- FALSE` en `00_configuracion.R`.

## Entorno

- R 4.2.1 en Windows. Paquetes principales: `dplyr`, `readr`, `tidyr`, `ggplot2`, `lubridate`, `stringr`, `scales`, `MASS`, `mgcv`, `gbm`, `ranger`.
- Python 3 en Google Colab para la construcción de los datos y el análisis espacial (`pandas`, `numpy`, `geopandas`, `rasterio`, `libpysal`, `esda`).

## Reproducibilidad

- Semilla de los modelos y del bootstrap: 2024. Semilla del análisis espacial: 12345. Ambas se fijan en el código.
- Partición temporal: entrenamiento 2015-2021, prueba 2022-2024.
- Regla anti-fuga: para predecir la semana *t* con horizonte *h* solo se usa información disponible hasta la semana *t − h*.
- Los resultados pueden variar en las últimas cifras decimales entre versiones de R y de los paquetes.

## Cita

> Zúñiga, N. (2026). *Accidente ofídico en Colombia, 2015-2024: clima, brotes y predicción* (Versión 1.0) [Software]. GitHub. https://github.com/nizuro/accidente-ofidico-colombia-clima
