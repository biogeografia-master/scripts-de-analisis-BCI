#' ---
#' title: "Medición de asociación. Modo R aplicado a mi familia asignada"
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
#' ## Modo R: matrices de dependencia entre variables (índice de correlación)
#' 
#' ### Modo R para datos cuantitativos de especies (abundancia)
#' 
#' En este caso, las variables usaré los valores de abundancias de especies como variables. Es decir, **compararé el grado de asociación entre especies, NO entre sitios**.
#' 
#' Aunque se podría usar el índice de correlación como métrica de la dependencia (tal como mostré en el script `aed_5_correlaciones_variables_ambientales.R`), tanto los doble-ceros (ausencias de una misma especie en dos lugares), como los *outliers* de abundancia, contribuirían a aumentar de manera ficticia el índice de correlación de Pearson.
#' 
#' Por tal razón, es recomendable aplicar la transformación *Chi* a la matriz de comunidad transpuesta, logrando así comparar especies no sitios. La transposición de la matriz de comunidad es necesaria para garantizar que lo que se compara son especies, no sitios (modo R). Explico el procedimiento a continuación, paso a paso.
#' 
#' Primero, sustituyo el caracter de espacio por un <enter> en los nombres de las especies (caracter \\n), para facilitar la lectura de los nombres en la diagonal del mapa de calor. Luego transpongo la matriz.
#' 
mi_fam_t <- mc_apcyn_melic_saptc %>% 
  rename_all(gsub, pattern = ' ', replacement = '\n') %>% 
  t()
#' 
#' Segundo, transformo la matriz transpuesta usando estandarización *Chi*.
#' 
mi_fam_t_chi <- decostand(mi_fam_t, "chi.square")
#' 
#' Tercero,  calculo la distancia euclídea.
#' 
mi_fam_t_chi_d <- dist(mi_fam_t_chi)
#' 
#' Finalmente, creo el "mapa de calor".
#' 
coldiss(mi_fam_t_chi_d, diag = TRUE)
#'
#' En el mapa de calor ordenado (el de la derecha), se identifica al menos un patrón de dependencia entre las especies listadas en la diagonal desde *Chrysophyllum cainito* hasta *Trichilia pallida* (todas relacionadas entre sí por medio de cuadros de color rosa). También se observan las especies que no parecen asociarse con otras, situadas en los extremos de la diagonal, y relacionadas con otras por medio de valores pequeños de distancia (cuadros azules), como *Rauvolfia littoralis* y *Pouteria fossicola* y *Cedrela odorata*.
#' 
#' ### Modo R para datos binarios (presencia/ausencia)
#' 
#' Arriba usé la distancia de Jaccard para evaluar asociación entre sitios. Dicha métrica también se puede usar para evaluar la distancia entre especies, usando como fuente la matriz de comunidad transpuesta convertida a binaria (presencia/ausencia)
#' 
mi_fam_t_jac <- vegdist(mi_fam_t, "jaccard", binary = TRUE)
coldiss(mi_fam_t_jac, diag = TRUE)
#'