#' ---
#' title: "Análisis espacial de datos ecológicos <br> Parte 1: Autocorrelación"
#' author: "JR"
#' date: "5 de diciembre, 2020"
#' output: github_document
#' ---
#'
knitr::opts_chunk$set(fig.width=12, fig.height=8)
#'
#' ## Preámbulo
#' 
#' ### Cargar paquetes
#' 
library(ape)
library(spdep)
library(ade4)
library(adegraphics)
library(adespatial)
library(vegan)
library(tidyverse)
library(sf)
source('biodata/funciones.R')
#' 
#' ### Cargar datos
#' 
load('biodata/Apocynaceae-Meliaceae-Sapotaceae.Rdata')
load('biodata/matriz_ambiental.Rdata')
mi_fam <- mc_apcyn_melic_saptc
bci_env_grid %>% tibble
bci_env_grid_sp <- bci_env_grid %>% as_Spatial
centroides <- bci_env_grid %>% st_centroid
bci_xy <- centroides %>% st_coordinates %>% as.data.frame
vecindad <- bci_env_grid_sp %>% poly2nb
#' 
#' ## Autocorrelación variable ambiental
#' 
mi_fam_hel <- decostand (mi_fam, "hellinger")
plot(bci_env_grid_sp)
plot(vecindad, coords = bci_xy, add=T, col = 'red')
var_ph <- bci_env_grid %>% st_drop_geometry %>% pull(pH)
ph_correl <- sp.correlogram(vecindad,
                            var_ph,
                            order = 9,
                            method = "I",
                            zero.policy = TRUE)
print(ph_correl, p.adj.method = 'holm')
plot(ph_correl)
#'
#' ## Autocorrelación densidad de individuos por especies
#' 
mi_fam %>% colSums %>% sort
auto_spp <- calcular_autocorrelacion_especies(
  mc = mi_fam,
  area_sitio = 100,
  orden = 9)
print(auto_spp, p.adj.method = 'holm')
dim_panel <- rev(n2mfrow(ncol(mi_fam)))
par(mfrow = dim_panel)
invisible(lapply(auto_spp, function(x) plot(x, main = x$var)))
#' 
#' Correlograma Mantel
#' 
mi_fam_sin_tendencia <- resid(
  lm(as.matrix(mi_fam_hel) ~ .,
     data = bci_xy))
mi_fam_sin_tendencia_d <- dist(mi_fam_sin_tendencia)
(mi_fam_correlograma <- mantel.correlog(
  mi_fam_sin_tendencia_d,
  XY = bci_xy,
  nperm = 999))
summary(mi_fam_correlograma)
plot(mi_fam_correlograma)
