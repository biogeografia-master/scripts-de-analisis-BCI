unir_censo_con_lista <- function(censo, lista) {
  library(dplyr)
  salida <- censo %>% inner_join(lista, by = 'sp')
  return(salida)
}

filtrar_vivos <- function(df){
  library(dplyr)
  salida <- df %>% filter(DFstatus=='alive' & dbh > 0)
  return(salida)
}

generar_campo_quad1ha <- function(df){
  library(dplyr)
  load('CTFSRPackage.rdata')
  salida <- df %>% mutate(quad1ha = gxgy.to.hectindex(gx, gy))
  return(salida)
}

filtrar_familia <- function(df, familia){
  library(dplyr)
  salida <- df %>% filter(Family %in% familia)
}

unir_lista_vivos <- function(censo, lista) {
  ccl <- unir_censo_con_lista(censo = censo, lista = lista)
  v <- filtrar_vivos(ccl)
  return(v)
}

unir_lista_vivos_q1ha <- function(censo, lista) {
  ulv <- unir_lista_vivos(censo = censo, lista = lista)
  q1ha <- generar_campo_quad1ha(ulv)
  return(q1ha)
}

unir_lista_vivos_q1ha_familia <- function(censo, lista, familia) {
  ulvq1ha <- unir_lista_vivos_q1ha(censo = censo, lista = lista)
  f <- filtrar_familia(ulvq1ha, familia = familia)
  return(f)
}

generar_matriz_comunidad_vivos_q1ha_familia <- function(
  censo, lista, familia) {
  library(dplyr)
  ulvq1hafam <- unir_lista_vivos_q1ha_familia(
    censo = censo, lista = lista, familia = familia)
  salida <- ulvq1hafam %>%
    mutate(n = 1) %>% 
    dplyr::select(quad1ha, Latin, n) %>%
    group_by(quad1ha, Latin) %>%
    summarise(N=sum(n)) %>%
    filter(!is.na(quad1ha)) %>%
    spread(Latin, N) %>%
    replace(., is.na(.), 0) %>%
    column_to_rownames('quad1ha')
  return(salida)
}

generar_matriz_comunidad_vivos_q1ha_global <- function(censo, lista) {
  library(dplyr)
  ulvq1ha <- unir_lista_vivos_q1ha(censo = censo, lista = lista)
  salida <- ulvq1ha %>%
    mutate(n = 1) %>% 
    dplyr::select(quad1ha, Latin, n) %>%
    group_by(quad1ha, Latin) %>%
    summarise(N=sum(n)) %>%
    filter(!is.na(quad1ha)) %>%
    spread(Latin, N) %>%
    replace(., is.na(.), 0) %>%
    column_to_rownames('quad1ha')
  return(salida)
}

generar_mc <- function(censo) {
  library(dplyr)
  salida <- censo %>%
    mutate(n = 1) %>% 
    dplyr::select(quad1ha, Latin, n) %>%
    group_by(quad1ha, Latin) %>%
    summarise(N=sum(n)) %>%
    filter(!is.na(quad1ha)) %>%
    spread(Latin, N) %>%
    replace(., is.na(.), 0) %>%
    column_to_rownames('quad1ha')
  return(salida)
}

generar_matriz_comunidad_vivos_q1ha_familia_como_taxon <- function(
  censo, lista) {
  library(dplyr)
  ulvq1ha <- unir_lista_vivos_q1ha(
    censo = censo, lista = lista)
  salida <- ulvq1ha %>%
    mutate(n = 1) %>% 
    dplyr::select(quad1ha, Family, n) %>%
    group_by(quad1ha, Family) %>%
    summarise(N=sum(n)) %>%
    filter(!is.na(quad1ha)) %>%
    spread(Family, N) %>%
    replace(., is.na(.), 0) %>%
    column_to_rownames('quad1ha')
  return(salida)
}

crear_grafico_mosaico_de_mc <- function(mc, tam_rotulo = 12){
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  para_gg <- mc %>%
    rownames_to_column('quad1ha') %>% mutate(quad1ha = as.numeric(quad1ha)) %>% 
    gather(especie, abundancia, -quad1ha) %>%
    mutate(quad1hafac = factor(quad1ha))
  p <- ggplot(
    para_gg,
    aes(x = quad1ha, y = especie, fill = abundancia, label = abundancia)
    ) +
    geom_tile(col='black') +
    scale_fill_gradient(
      trans = 'log1p',
      low = "white",
      high = "red") +
    geom_text(size = tam_rotulo/4) +
    scale_x_continuous(#For duplicate axis
      breaks = 1:length(levels(para_gg$quad1hafac)),
      labels = levels(para_gg$quad1hafac),
      sec.axis = dup_axis()
    ) +
    # scale_x_discrete(position = 'bottom') +
    theme_bw() +
    coord_equal() +
    theme(
      legend.position="none",
      text = element_text(size = tam_rotulo),
      panel.background = element_rect(fill = 'white', colour = 'black'),
      panel.grid.major = element_line(colour = "grey", linetype = "dashed", size = 0.25),
      axis.text = element_text(colour="black")
    )
    return(p)
}

coldiss <- function(D, nc = 4, byrank = TRUE, diag = FALSE)
{
  require(gclus)

  if (max(D)>1) D <- D/max(D)
  
  if (byrank) {
    spe.color <- dmat.color(1-D, cm.colors(nc))
  }
  else {
    spe.color <- dmat.color(1-D, byrank=FALSE, cm.colors(nc))
  }
  
  spe.o <- order.single(1-D)
  speo.color <- spe.color[spe.o, spe.o]
  
  op <- par(mfrow=c(1,2), pty="s")
  
  if (diag) {
    plotcolors(spe.color, rlabels=attributes(D)$Labels, 
               main="Dissimilarity Matrix", 
               dlabels=attributes(D)$Labels)
    plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o], 
               main="Ordered Dissimilarity Matrix", 
               dlabels=attributes(D)$Labels[spe.o])
  }
  else {
    plotcolors(spe.color, rlabels=attributes(D)$Labels, 
               main="Dissimilarity Matrix")
    plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o], 
               main="Ordered Dissimilarity Matrix")
  }
  
  par(op)
  #example:
  #coldiss(matriz_dis, diag = TRUE)
}

guardar_rdata <- function(sufijo){
  load('abreviaturas_de_familias.rdata')
  save(
    list = grep(sufijo, ls(envir = .GlobalEnv), value = T),
    file = paste0(famabr[famabr$abr==sufijo, 'fam'], '.Rdata')
  )
}

organizar_matriz_distancia <- function(matriz_distancia, func_dist){
  library(dplyr)
  library(broom)
  matriz_distancia %>% 
    tidy() %>% 
    arrange(item2) %>% 
    mutate(
      combinacion = paste0('D(', paste(item2, item1, sep = ', '), ')'),
      `Funcion de distancia` = func_dist) %>% 
    select(-contains('item')) %>% 
    spread(combinacion, distance)
}

ezonalobj <- function(objraster = NULL, nombre = '', objgeometrias = NULL, export = T, nombreexport = '', cuali = F){
  #Ejemplo cuantitativo: ezonal(objraster = 'vrm_90M_n00w090/vrm_mosaico.tif', nombre = 'vrm', objgeometrias = 'divisionRD.gpkg', capa = 'MUNCenso2010', export = T, nombreexport = 'divisionRD_vrm', cuali = F)
  #Ejemplo cualitativo: ezonal(objraster = 'geom_90M_n00w090/geom_mosaico.tif', nombre = 'geomorfonos', objgeometrias = 'divisionRD.gpkg', capa = 'MUNCenso2010', export = T, nombreexport = 'divisionRD_geom', cuali = T)
  require(raster)
  require(sf)
  require(dplyr)
  require(tidyr)
  r <- objraster
  names(r) <- nombre
  geoms <- objgeometrias
  df <- raster::extract(x = r, y = as(geoms, 'Spatial'), df = T)
  if(cuali){
    dfresumen <- df %>%
      count(.[[1]], .[[2]]) %>%
      group_by(.[[1]]) %>%
      mutate(pct=n/sum(n)*100) %>%
      dplyr::select(-n) %>%
      spread(`.[[2]]`, pct) %>%
      ungroup() %>%
      dplyr::select(-`.[[1]]`) %>%
      rename_all(function(x) paste(nombre, x, sep = '_'))
  } else {
    dfresumen <- t(sapply(unique(df[,'ID']), function(x) c(
      n = length(na.omit(df[df[,'ID']==x,nombre])),
      min = min(df[df[,'ID']==x,nombre], na.rm = T),
      cuartil_ = quantile(df[df[,'ID']==x,nombre], 1/4, na.rm = T),
      media = mean(df[df[,'ID']==x,nombre], na.rm = T),
      mediana = median(df[df[,'ID']==x,nombre], na.rm = T),
      cuartil_ = quantile(df[df[,'ID']==x,nombre], 3/4, na.rm = T),
      max = max(df[df[,'ID']==x,nombre], na.rm = T),
      desv = sd(df[df[,'ID']==x,nombre], na.rm = T))))
    colnames(dfresumen) <- paste0(nombre, '_', colnames(dfresumen))
    colnames(dfresumen) <- gsub('\\.', '', colnames(dfresumen))
  }
  # return(dfresumen)
  geomsout <- dplyr::bind_cols(geoms, as.data.frame(dfresumen))
  if (export) {
    st_write(geomsout, dsn = paste0('salidas_ezonal/', nombreexport, '.gpkg'), driver = 'GPKG')
    saveRDS(geomsout, file = paste0('salidas_ezonal/', nombreexport,'.RDS'))
  }
  return(geomsout)
}

coldissgg <- function(dist, ordered = T, nc = 100, fsz = 4) {
  #dist = matriz de distancias
  #ordered = ordered by distance value
  #nc = number of colors
  require(reshape2)
  require(tidyr)
  require(dplyr)
  require(gclus)
  require(RColorBrewer)
  dist.g <- melt(as.matrix(dist))
  dist.g[,c('Var1','Var2')] <- lapply(dist.g[,c('Var1','Var2')], factor)
  dist.g$Var1 <- factor(dist.g$Var1, levels=sort(levels(dist.g$Var1)))
  dist.g$Var2 <- factor(dist.g$Var2, levels=sort(levels(dist.g$Var2)))
  dist.g[dist.g$Var1==dist.g$Var2,'value'] <- NA
  dist.g$type <- 'Dissimilarity matrix'
  dist.g.o <- dist.g
  dist.g.o$Var1 <- factor(dist.g.o$Var1, levels=levels(dist.g.o$Var1)[order.single(1-dist)])
  dist.g.o$Var2 <- factor(dist.g.o$Var2, levels=levels(dist.g.o$Var2)[order.single(1-dist)])
  dist.g.o$type <- 'Ordered dissimilarity matrix'
  mypalette <- colorRampPalette(rev(brewer.pal(4, "YlOrRd")), space="Lab")#Borrowed from Heatmap.R with spectral palette: https://gist.github.com/dsparks/3710171
  gg1 <- ggplot(dist.g, aes(Var1, Var2)) +
    geom_tile(aes(fill=value), colour = "white") +
    scale_fill_gradientn(colours = mypalette(nc), na.value = 'white') +
    geom_text(aes(label=round(value,2)), size = fsz) +
    labs(title='Dissimilarity matrix') +
    theme(
      text = element_text(size = 16),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position = 'none',
      plot.title = element_text(size=18, hjust = 0.5)
      ) +
    coord_equal()
  gg2 <- ggplot(dist.g.o, aes(Var1, Var2)) +
    geom_tile(aes(fill=value), colour = "white") +
    scale_fill_gradientn(colours = mypalette(nc), na.value = 'white') +
    geom_text(aes(label=round(value,2)), size = fsz) +
    labs(title='Ordered dissimilarity matrix') +
    theme(
      text = element_text(size = 16),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position = 'none',
      plot.title = element_text(size=18, hjust = 0.5)
    ) +
    coord_equal()
  if(ordered) print(gg2) else print(gg1)
}
