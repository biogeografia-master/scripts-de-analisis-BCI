#' ---
#' title: "Medición de asociación. Modo Q aplicado a mi familia asignada"
#' author: "JR"
#' date: "3 de noviembre, 2020"
#' output: github_document
#' ---

knitr::opts_chunk$set(fig.width=12, fig.height=8)

#' ## Preámbulo

#' ### Cargar paquetes
library(ez)
library(psych)
library(vegan)
library(adespatial)
library(broom)
library(tidyverse)
library(sf)
source('biodata/funciones.R')

#' ### Cargar datos
#' 
load('biodata/matriz_ambiental.Rdata')
load('biodata/Apocynaceae-Meliaceae-Sapotaceae.Rdata')
#' 
#' ## Modo Q: matrices de disimilaridad entre objetos
#' 
#' ### Modo Q para datos cuantitativos de especies (abundancia). Datos de mi familia asignada
#' 
#' Aplicado a BCI y mi familia asignada (en forma de matriz de distncia), utilizando la transformación *Hellinger*:
#' 
(mi_fam_d_hel <- dist.ldc(mc_apcyn_melic_saptc, "hellinger", silent = T))
#' 
#' Para interpretar esta matriz, es necesario representarla gráficamente. Color fucsia (magenta, rosa) significa "corta distancia=muy similares", y cian (celeste) significa "gran distancia=poco similares:
#' 
coldiss(mi_fam_d_hel, diag = T)
#' 
#' Puedes guardar el gráfico usando el botón `Export` de la pestaña `Plots`
#' 
#' Una forma alterna de guardar el gráfico es mediante funciones de R. La calidad de gráficos exigida en revistas, suele requerir usar dichas funciones específicas, porque permiten más control. Abajo uso la función `png` para "abrir un dispositivo gráfico", luego imprimo el gráfico que deseo guardar y finalmente cierro el dispositivo mediante `dev.off` Por ejemplo:
#' 
png(
  filename = 'matriz_disimilaridad_hellinger.png',
  width = 2400, height = 1200, pointsize = 32
)
coldiss(mi_fam_d_hel, diag = T)
dev.off()
#' 
#' MUY IMPORTANTE. La última función, `dev.off()`, es necesaria para cerrar el dispositivo. Si no la ejecutas, no se generarán gráficos en el dispositivo estándar (e.g. pestaña `Plots`) 
#' 
#' ### Modo Q para datos binarios (presencia/ausencia)
#' 
#' Frecuentemente, sólo dispones de datos de presencia/ausencia. En esos casos, existe un conjunto de herramientas con las que analizar posibles patrones de asociación. A continuación muestro cómo calcular la distancia de Jaccard.
#' 
mi_fam_jac <- vegdist(mc_apcyn_melic_saptc, method = 'jac', binary = T)
mi_fam_jac
#' 
#' El argumento `binary=T` "ordena" que se realice primero `decostand(mc_apcyn_melic_saptc, method = 'pa')`, y se convierte a la matriz de comunidad en una de presencia/ausencia.
#' 
#' En esta matriz de disimilaridad, al igual que en la anterior, un valor pequeño significa que los sitios comparados son muy parecidos.
coldiss(mi_fam_jac, diag = T)
#' 
#' A continuación, usando la distancia Sorensen o Bray-Curtis:
#' 
mi_fam_sor <- vegdist(mc_apcyn_melic_saptc, method = 'bray', binary = T)
coldiss(mi_fam_sor, diag = T)
#' 
#' ### Modo Q para datos cuantitativos, excluyendo abundancia de especies (variables ambientales)
#' 
#' En este ejemplo, usaré sólo variables de suelo, todas cuantitativas, puedes combinar con otras variables que hayas detectado como relevantes en el análisis de correlación
#' 
env_suelo_punt_z <- bci_env_grid %>%
  st_drop_geometry() %>% 
  dplyr::select(matches('^[A-T,Z]|^pH$', ignore.case = F)) %>% 
  scale()
env_suelo_punt_z_d <- dist(env_suelo_punt_z)
coldiss(env_suelo_punt_z_d, diag = T)
#'
#' ### Modo Q para datos cualitativos y cuantitativos (mixtos), excluyendo abundancia de especies (variables ambientales)
#' En este ejemplo, usaré las siguientes variables mixtas (funciona igualmente para datos cualitativos solamente):
#' 
#' - `hetereogeneidad_ambiental`. Índice cuantitativo que informa sobre la heterogeneidad basada en múltiples variables topográficas y composición de especies.
#' 
#' 
#' - `habitat`. Tipo de hábitat. Asume los siguientes valores posibles: *OldHigh*, *OldLow* y *OldSlope* (bosque viejo en relieve alto, en vertientes y relieve bajo, respectivamente), *Swamp* (bosque en área encharcable) *Young* (bosque joven).
#' 
#' - `quebrada`. Informa sobre si hay o no quebrada. Los valores posibles son *Yes* o *No*.
#' 
env_mix <- bci_env_grid %>%
  st_drop_geometry() %>%
  dplyr::select(heterogeneidad_ambiental, habitat, quebrada)
env_mix_d <- daisy(x = env_mix, metric = 'gower')
env_mix_d %>% coldiss(diag = T)
#'