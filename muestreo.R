#' ---
#' title: "Muestreo"
#' author: "JR"
#' date: "3 de noviembre, 2020"
#' output: github_document
#' ---

knitr::opts_chunk$set(fig.width=12, fig.height=8)

#' ## Cargar paquetes
library(raster)
library(mapview)
library(RColorBrewer)
library(sf)

#' ## Lectura de datos fuente, generación de ráster simplicado
#' Tree cover, from here: http://earthenginepartners.appspot.com/google.com/science-2013-global-forest
# r <- raster('treecover2000_crop.tif')
# r2 <- raster::aggregate(r, 20)
# r2
# writeRaster(r2, 'treecover2000_remuestreado.tif')

#' ## Lectura de datos procesados (ráster simplificado)
r2 <- raster(x = 'treecover2000_remuestreado.tif')

#' ## Mapa (exploratorio)
verdes <- colorRampPalette(brewer.pal(8, "Greens"))
mapview(r2, col.regions = verdes, alpha = 0.5)

#' ## Reclasificación del ráster (dos clases, 0-25%=clase 1, 26-100%=clase 2)
m <- c(-1, 25, 1,  25, 100, 2)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- reclassify(r2, rclmat)
rc

#' ## Desplegar en mapview
mapview(rc, method = 'ngb')

#' ## Conversión a polígonos
# poligonos <- rasterToPolygons(rc)
#' Guardar polígonos
# saveRDS(poligonos, 'poligonos_raster_cobertura.RDS')
poligonos <- readRDS('poligonos_raster_cobertura.RDS')

#' ## Tomar muestra estratificada
set.seed(10)
(muestra <- spsample(poligonos, n = 100, type = 'stratified'))
# saveRDS(muestra, 'muestra.RDS')

#' Asignarle atributos a cada punto de muestra
muestra_con_atr <- over(muestra, poligonos)
muestra$clase_de_cobertura <- muestra_con_atr

#' Convertir a sf (necesario para representar con mapview)
muestra_sf <- st_as_sf(muestra)

#' Representar
mapview(rc, method = 'ngb', alpha = 0.5) + mapview(muestra_sf)

#' Guardar muestra_sf
saveRDS(muestra_sf, 'muestra_sf.RDS')
