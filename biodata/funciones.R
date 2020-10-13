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
