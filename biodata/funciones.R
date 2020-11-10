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


# panelutils.R 
#
# License: GPL-2 
# Author: Francois Gillet, 23 August 2012
#
## Put Pearson, Spearman or Kendall correlations on the upper panel
panel.cor <- function(x, y, method="pearson", digits=3, cex.cor=1.2, no.col=FALSE)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, method=method)
  ra <- cor.test(x, y, method=method)$p.value
  txt <- round(r, digits)
  prefix <- ""
  if(ra <= 0.1) prefix <- "."
  if(ra <= 0.05) prefix <- "*"
  if(ra <= 0.01) prefix <- "**"
  if(ra <= 0.001) prefix <- "***"
  if(no.col)
  {
    color <- 1
    if(r < 0) { if(ra <= 0.001) sig <- 4 else sig <- 3 }
    else { if(ra <= 0.001) sig <- 2 else sig <- 1 }
  }
  else
  {
    sig <- 1
    if(ra <= 0.001) sig <- 2
    color <- 2
    if(r < 0) color <- 4
  }
  txt <- paste(txt, prefix, sep="\n")
  text(0.5, 0.5, txt, cex = cex.cor, font=sig, col=color)
}


## Put histograms on the diagonal
panel.hist <- function(x, no.col=FALSE, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  his <- hist(x, plot=FALSE)
  breaks <- his$breaks; nB <- length(breaks)
  y <- his$counts
  y <- y/max(y)
  if(no.col) rect(breaks[-nB], 0, breaks[-1], y, col="gray", ...)
  else rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}


## Add black lowess curves to scatter plots
panel.smoothb <- function (x, y, col=par("col"), bg=NA, pch=par("pch"), 
                           cex=1, col.smooth="black", span=2/3, iter=3, ...) 
{
  points(x, y, pch=pch, col=col, bg=bg, cex=cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f=span, iter=iter), col=col.smooth, ...)
}
#Usage:
#pairs(num.mat, lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=panel.hist)
#pairs(num.mat, lower.panel=panel.smooth, upper.panel=panel.cor, method="kendall")

ezCorM <- function (data, r_size_lims = c(10, 30), point_alpha = 0.5, density_height = 1, 
                    density_adjust = 1, density_colour = "white", label_size = 10, 
                    label_colour = "black", label_alpha = 0.5, lm_colour = "red", 
                    ci_colour = "green", ci_alpha = 0.5, test_alpha = 0.05, test_correction = "none", method = 'pearson') 
{
  library(plyr)
  library(reshape2)
  library(ggplot2)
  if (inherits(data, "tbl_df")) {
    data <- as.data.frame(data)
  }
  ntests = ((((ncol(data) - 1)^2) - (ncol(data) - 1))/2)
  if (test_correction[1] == "bonferroni") {
    test_alpha = test_alpha/ntests
  }
  else {
    if (test_correction[1] == "sidak") {
      test_alpha = 1 - (1 - test_alpha)^(1/ntests)
    }
  }
  for (i in 1:length(data)) {
    data[, i] = (data[, i] - mean(data[, i], na.rm = T))/sd(data[, 
                                                                 i], na.rm = T)
  }
  z = data.frame()
  z_cor = data.frame()
  i = 1
  j = i
  while (i <= length(data)) {
    if (j > length(data)) {
      i = i + 1
      j = i
    }
    else {
      x = data[, i]
      y = data[, j]
      toss = is.na(x) | is.na(y)
      x = x[!toss]
      y = y[!toss]
      temp = as.data.frame(cbind(x, y))
      temp = cbind(temp, names(data)[i], names(data)[j])
      z = rbind(z, temp)
      this_cor = round(cor(x, y, method = method), 2)
      this_cor.test = cor.test(x, y, method = method, exact = F)
      this_col = ifelse(this_cor.test$p.value < test_alpha, 
                        "a", "b")
      this_size = (this_cor)^2
      cor_text = ifelse(this_cor == 0, "0", ifelse(this_cor == 
                                                     1, "1", ifelse(this_cor == -1, "-1", ifelse(this_cor > 
                                                                                                   0, substr(format(c(this_cor, 0.123456789), digits = 2)[1], 
                                                                                                             2, 4), paste("-", substr(format(c(this_cor, 0.123456789), 
                                                                                                                                             digits = 2)[1], 3, 5), sep = "")))))
      b = as.data.frame(cor_text)
      b = cbind(b, this_col, this_size, names(data)[j], 
                names(data)[i])
      z_cor = rbind(z_cor, b)
      j = j + 1
    }
  }
  names(z) = c("x", "y", "x_lab", "y_lab")
  z = z[z$x_lab != z$y_lab, ]
  names(z_cor) = c("cor", "p", "rsq", "x_lab", "y_lab")
  z_cor = z_cor[z_cor$x_lab != z_cor$y_lab, ]
  diag = melt(data, measure.vars = names(data))
  names(diag)[1] = "x_lab"
  diag$y_lab = diag$x_lab
  dens = ddply(diag, .variables = .(x_lab, y_lab), function(x) {
    d = density(x$value[!is.na(x$value)], adjust = density_adjust)
    d = data.frame(x = d$x, y = d$y)
    d$ymax = d$y * (max(abs(c(z$x, z$y))) * 2 * density_height)/max(d$y) - 
      max(abs(c(z$x, z$y))) * density_height
    d$ymin = -max(abs(c(z$x, z$y))) * density_height
    return(d)
  })
  labels = ddply(diag, .variables = .(x_lab, y_lab), function(x) {
    to_return = data.frame(x = 0, y = 0, label = x$x_lab[1])
    return(to_return)
  })
  points_layer = layer(geom = "point", params = list(alpha = point_alpha, 
                                                     na.rm = TRUE), stat = "identity", position = "identity", 
                       data = z, mapping = aes_string(x = "x", y = "y"))
  lm_line_layer = layer(geom = "line", stat = "smooth", position = "identity", 
                        params = list(colour = lm_colour, method = "lm", na.rm = TRUE), 
                        data = z, mapping = aes_string(x = "x", y = "y"))
  lm_ribbon_layer = layer(geom = "ribbon", stat = "smooth", 
                          position = "identity", params = list(fill = ci_colour, 
                                                               alpha = ci_alpha, method = "lm", na.rm = TRUE), data = z, 
                          mapping = aes_string(x = "x", y = "y"))
  cor_text_layer = layer(geom = "text", stat = "identity", 
                         position = "identity", data = z_cor, mapping = aes_string(label = "cor", 
                                                                                   size = "rsq", colour = "p"), params = list(x = 0, 
                                                                                                                              y = 0, na.rm = TRUE))
  dens_layer = layer(geom = "ribbon", stat = "identity", position = "identity", 
                     params = list(colour = "transparent", fill = "white", 
                                   na.rm = TRUE), data = dens, mapping = aes_string(x = "x", 
                                                                                    ymax = "ymax", ymin = "ymin"))
  label_layer = layer(geom = "text", stat = "identity", position = "identity", 
                      params = list(colour = label_colour, size = label_size, 
                                    alpha = label_alpha, na.rm = TRUE), data = labels, 
                      mapping = aes_string(x = "x", y = "y", label = "label"))
  y_lab = NULL
  x_lab = NULL
  f = facet_grid(y_lab ~ x_lab)
  o = theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
            axis.ticks = element_blank(), axis.text.y = element_blank(), 
            axis.text.x = element_blank(), axis.title.y = element_blank(), 
            axis.title.x = element_blank(), legend.position = "none", 
            strip.background = element_blank(), strip.text.x = element_blank(), 
            strip.text.y = element_blank())
  x_scale = scale_x_continuous(limits = c(-1 * max(abs(dens$x)), 
                                          max(abs(dens$x))))
  size_scale = scale_size(limits = c(0, 1), range = r_size_lims)
  return(ggplot(z_cor) + points_layer + lm_ribbon_layer + lm_line_layer + 
           dens_layer + label_layer + cor_text_layer + f + o + x_scale + 
           size_scale)
}