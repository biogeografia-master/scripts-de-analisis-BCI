#' ---
#' title: "Cómo integrar varias imágenes (e.g. gráficos, fotos, mapas) en una sola (aka 'panel')"
#' author: "JR"
#' date: "8 de enero, 2021"
#' output: github_document
#' ---

knitr::opts_chunk$set(fig.width=12, fig.height=8)

#'  ## Cargar paquetes
library(cowplot)
library(ggplot2)
library(magick) #En caso de error, "sudo apt install libmagick++-dev"
#' 
#' ## Imágenes
#' 
pngs <- list.files(pattern='*.png')
sample(pngs, 3)
#' 
#' ## Creación repetitiva (sorry, "do not repeat yourself")
#' 
p1 <- ggdraw() + draw_image(pngs[1], scale = 0.9)
p2 <- ggdraw() + draw_image(pngs[2], scale = 0.9)
p3 <- ggdraw() + draw_image(pngs[3], scale = 0.9)
#' 
#' ## Generar JPG
#' 
jpeg('ejemplo_de_panel_varios_plots.jpg', width = 1920, height = 540, res = 150)
plot_grid(p1, p2, p3, labels = c('(a)','(b)','(c)'), nrow = 1)
dev.off()
#'
#' ## Usando función `crear_panel`
#' 
#' Primero, carga la función:
#' 
#' Desde un clon del repo `scripts-de-analisis-BCI`
#' 
source('biodata/funciones.R')
#'
#' Desde cualquier lugar
#' 
devtools::source_url('https://raw.githubusercontent.com/biogeografia-master/scripts-de-analisis-BCI/master/biodata/funciones.R')
#' Ejecútala
crear_panel('mapa_cuadros_pendiente.png', 'mapa_cuadros_abun_global.png')
#'
#' Si prefieres guardar en archivo de salida
#' 
jpeg('ejemplo_de_panel_varios_plots.jpg', width = 1080, height = 540, res = 150)
crear_panel('mapa_cuadros_pendiente.png', 'mapa_cuadros_abun_global.png')
dev.off()
