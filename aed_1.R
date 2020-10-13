#' ---
#' title: "Análisis exploratorio de datos, parte 1"
#' author: "JR"
#' date: "13 de octubre, 2020"
#' output: github_document
#' ---

#' ### Área de cargar paquetes
library(vegan)
library(tidyverse)
library(sf)
source('biodata/funciones.R')

#' ### Área de cargar datos
#' Censo (el objeto se carga con prefijo "censo") y matriz de comunidad (prefijo "mc")
load('biodata/Apocynaceae-Meliaceae-Sapotaceae.Rdata')
load('biodata/matriz_ambiental.Rdata') #Matriz ambiental, se carga como "bci_env_grid"

#' ### Imprimir datos en pantalla (impresiones parciales con head)
head(censo_apcyn_melic_saptc)
head(mc_apcyn_melic_saptc)
bci_env_grid # No necesita imprimirse parcialmente

#' ### También podemos usar
#' Requiere que se haya cargado ya la colección tidyverse
censo_apcyn_melic_saptc %>% tibble
mc_apcyn_melic_saptc %>% tibble

#' ### Lista de especies
sort(colnames(mc_apcyn_melic_saptc))

#' ### Número de sitios, tanto en matriz de comunidad como en ambiental
#' Verifica que coinciden
nrow(mc_apcyn_melic_saptc) #En la matriz de comunidad
nrow(bci_env_grid) #En la matriz ambiental

#' ### Riqueza numérica de especies (usando matriz de comunidad) por quadrat
#' Nota: cargar paquete vegan arriba, en el área de paquetes
specnumber(mc_apcyn_melic_saptc)
sort(specnumber(mc_apcyn_melic_saptc)) # Ordenados ascendentemente
summary(specnumber(mc_apcyn_melic_saptc)) # Resumen estadístico

#' ### Abundancia de especies por quadrat
sort(rowSums(mc_apcyn_melic_saptc))
summary(rowSums(mc_apcyn_melic_saptc)) # Resumen estadístico

#' ### Abundancia por especie
sort(colSums(mc_apcyn_melic_saptc))
summary(colSums(mc_apcyn_melic_saptc)) # Resumen estadístico

#' ### Riqueza numérica de toda la "comunidad"
specnumber(colSums(mc_apcyn_melic_saptc))

#' ### Abundancia de toda la comunidad
sum(colSums(mc_apcyn_melic_saptc))

#' ### Una tabla para el manuscrito, es necesario asignarle nombre
#' Para esto, usaré la colección "tidyverse"
abun_sp <- censo_apcyn_melic_saptc %>%
  group_by(Latin) %>% 
  count() %>% 
  arrange(desc(n))
abun_sp

#' ### Un gráfico para el manuscrito
#' Gráfico de mosaicos de la abundancia por especie por cuadros
abun_sp_q <- crear_grafico_mosaico_de_mc(mc_apcyn_melic_saptc, tam_rotulo = 6)
abun_sp_q