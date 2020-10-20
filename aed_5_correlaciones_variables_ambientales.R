#' ---
#' title: "Análisis exploratorio de datos. Correlaciones entre variables ambientales"
#' author: "JR"
#' date: "18 de octubre, 2020"
#' output: github_document
#' ---

#' ### Cargar paquetes
library(tidyverse)
library(sf)
library(ez)
library(psych)

#' ### Cargar datos
load('biodata/matriz_ambiental.Rdata')

#' ### Una correlación simple
cor(bci_env_grid$pendiente_media, bci_env_grid$geomorf_vertiente_pct)
plot(bci_env_grid$pendiente_media, bci_env_grid$geomorf_vertiente_pct)
cor.test(bci_env_grid$pendiente_media, bci_env_grid$geomorf_vertiente_pct)

#' ### Generar objeto de columnas numéricas
#' #### Alternativa 1
env_num <- st_drop_geometry(bci_env_grid[sapply(bci_env_grid, is.numeric)])
#' #### Alternativa 2
env_num <- bci_env_grid %>% select_if(is.numeric) %>% st_drop_geometry

#' ### Evaluando correlación
pairs(env_num)
cor(env_num)
pairs.panels(env_num)

#' ### Panel de correlaciones
#' #### Todas las variables (se empastan)
p_cor_todos <- bci_env_grid %>%
  select_if(is.numeric) %>%
  st_drop_geometry() %>%
  dplyr::select(-id) %>%
  ezCor(r_size_lims = c(4,8), label_size = 4)
p_cor_todos
#' #### Sin "geomórfonos"
p_cor_sin_geomorf <- bci_env_grid %>%
  select_if(is.numeric) %>%
  st_drop_geometry() %>%
  dplyr::select(-id, -contains('geomorf')) %>%
  ezCor(r_size_lims = c(4,8), label_size = 4)
p_cor_sin_geomorf
#' #### Riqueza y abundancia con "geomórfonos"
p_cor_geomorf_ar <- bci_env_grid %>%
  select_if(is.numeric) %>%
  st_drop_geometry() %>%
  dplyr::select(matches('geomorf|abundancia|riqueza')) %>%
  ezCor(r_size_lims = c(4,8), label_size = 4)
p_cor_geomorf_ar
