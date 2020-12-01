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

calcular_anchuras_siluetas <- function(mc_orig, distancias, cluster) {
  siluetas <- numeric(nrow(mc_orig)) # Vector vacío, que acogera las siluetas
  i <- NULL; sil_i <- NULL
  for (i in 2:(nrow(mc_orig) - 1)) {
    sil_i <- silhouette(cutree(cluster, k = i), distancias)
    siluetas[i] <- summary(sil_i)$avg.width
  }
  siluetas
  n_grupos_optimo <- which.max(siluetas)
  return(
    list(
      anchuras_siluetas = siluetas,
      n_grupos_optimo = n_grupos_optimo
    )
  )
}

'cleanplot.pca' <- 
  function(res.pca, ax1=1, ax2=2, scaling=2, plot.sites=TRUE, 
           plot.spe=TRUE, label.sites=TRUE, label.spe=TRUE, cex.char1=0.7,
           pos.sites=2, pos.spe=4, mult.spe=1, select.spe=NULL, 
           mar.percent=0.1, optimum=TRUE, move.origin=c(0,0), silent=TRUE)
    # FUENTE: https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/cleanplot.pca.R
    # A function to draw a triplot (scaling 1 or scaling 2) from an object 
    # of class "rda" containing RDA result from vegan's rda() function.
    #
    # ARGUMENTS
    #
    # ##### General parameters
    # res.pca          An rda{vegan} object.
    # ax1, ax2         Canonical axes to be drawn as abscissa and ordinate. Defaults: 1 and 2.
    # site.sc          Can be set to "lc" (linear constraints or model scores, default) 
    #                  or "wa" (weighted averages, default in vegan).
# scaling          Scaling type: only 1 or 2 are supported. Default: 2.
#
# ##### Items to be plotted
# plot.sites       If TRUE, the sites will be plotted as small circles.
# plot.spe         If TRUE, the species (or other response variables) will be plotted.
# label.sites      If TRUE, labels are added to the site symbols.
# label.spe        If TRUE, labels are added to the species arrows.
# label.env        If TRUE, labels are added to the environmental variable arrows.
# label.centr      If TRUE, labels are added to the centroids of factor levels.
# cex.char1        Character size (for sites and response variables).
#
# ##### Label positions
# ## Positions: 1=below the point, 2=left, 3=above, 4=right. Default: 4.
# ## Note - Argument pos=NULL centres the label on the position of the object (site point,  
# ## species or environmental variable arrow, centroid) when the object is not drawn.
# pos.sites        Position of site labels. 1 to 4, as above. Default: 2.
# pos.spe          Position of species labels. 1 to 4, as above. Default: 4.
#
# ##### Multipliers, selection of species to be plotted
# mult.spe         Multiplier for length of the species arrows. Default: 1.
# select.spe       Vector containing a selection of the species numbers to be drawn in 
#                  the biplot, e.g. c(1,2,5,8,12). Draw all species if select.spe=NULL 
#                  (default value). The species that are well represented in the RDA plot 
#                  can be identified using goodness(RDA.output.object,display="species")
#
# ##### Position of the plot in frame, margins
# mar.percent      Factor to expand plot size to accomodate all items and labels. Positive 
#                  values increase the margins around the plot, negative values reduce 
#                  them.
# optimum          If TRUE, the longest species and environmental arrows are stretched to 
#                  a length equal to the distance to the origin of the site farthest from 
#                  the origin in the plot of (ax1,ax2). This is an optimal combined 
#                  representation of the three elements. The lengths of the species and 
#                  environmental arrows can be further modified using the arguments 
#                  mult.spe and mult.arrow.
# move.origin      Move plot origin right-left and up-down. Default: move.origin=c(0,0).
#                  Ex. move.origin=c(-1,0.5) moves origin by 1 unit left and 0.5 unit up.
#
# ##### Varia
# silent           If FALSE, intermediate computation steps are printed. Default: TRUE.
#
# # Example 1 - Table 11.3 of Legendre & Legendre (2012, p. 644), first 6 species only
#
# Y.mat = matrix(c(1,0,0,11,11,9,9,7,7,5,0,0,1,4,5,6,7,8,9,10,0,0,0,0,17,0,13,0,10,0,0, 
# 0,0,0,7,0,10,0,13,0,0,0,0,8,0,6,0,4,0,2,0,0,0,1,0,2,0,3,0,4),10,6)
# Depth = 1:10
# Sub. = as.factor(c(rep(1,3),4,2,4,2,4,2,4))
# env = cbind(data.frame(Depth),data.frame(Sub.))
# 
# rda.out = rda(Y.mat~ .,env)
# 
# # Scaling=1
# par(mfrow=c(1,2))
# triplot.rda(rda.out, scaling=1, mar.percent=0)
# triplot.rda(rda.out, scaling=1, move.origin=c(5,-5), mar.percent=-0.1)
#
# # Scaling=2
# par(mfrow=c(1,2))
# triplot.rda(rda.out, scaling=2, mar.percent=0.15, silent=FALSE)
# triplot.rda(rda.out, scaling=2, move.origin=c(0.4,-0.25), mar.percent=0.05,silent=FALSE)
#
# # Example 2 - Dune data
# 
# library(vegan)
# data(dune)
# data(dune.env)
# 
# rda.dune = rda(dune ~ .,dune.env)
# 
# tmp = goodness(rda.dune)
# ( sp.sel = which(tmp[,2] >= 0.4) )
#
# Scaling=2
# par(mfrow=c(1,2))
# triplot.rda(rda.dune, mar.percent=0)
# triplot.rda(rda.dune, select.spe=sp.sel, move.origin=c(-0.3,0), mar.percent=0.1)
#
# #####
#
# License: GPL-2 
# Authors: Francois Gillet, Daniel Borcard & Pierre Legendre, 2016
{
  ### Internal functions
  #
  'stretch' <- 
    function(sites, mat, ax1, ax2, n, silent=silent) {
      # Compute stretching factor for the species arrows
      # First, compute the longest distance to centroid for the sites
      tmp1 <- rbind(c(0,0), sites[,c(ax1,ax2)])
      D <- dist(tmp1)
      target <- max(D[1:n])
      # Then, compute the longest distance to centroid for the species arrows
      if(class(mat)=="matrix") {
        p <- nrow(mat)   # Number of species to be drawn
        tmp2 <- rbind(c(0,0), mat[,c(ax1,ax2)])
        D <- dist(tmp2)
        longest <- max(D[1:p])
      } else { tmp2 <- rbind(c(0,0), mat[c(ax1,ax2)]) 
      longest <- dist(tmp2)
      # print(tmp2)
      }  # If a single row left in 'mat'
      #
      if(!silent) cat("target =",target," longest =",longest," fact =",target/longest,"\n")
      fact <- target/longest
    }
  #
  'larger.plot' <- 
    function(sit.sc, spe.sc, percent, move.origin, ax1, ax2) {
      # Internal function to expand plot limits (adapted from code by Pierre Legendre)
      mat <- rbind(sit.sc, spe.sc)
      range.mat <- apply(mat, 2, range)
      rownames(range.mat) <- c("Min","Max")
      z <- apply(range.mat, 2, function(x) x[2]-x[1])
      range.mat[1,] <- range.mat[1,]-z*percent
      range.mat[2,] <- range.mat[2,]+z*percent
      if(move.origin[1] != 0) range.mat[,ax1] <- range.mat[,ax1] - move.origin[1]
      if(move.origin[2] != 0) range.mat[,ax2] <- range.mat[,ax2] - move.origin[2]
      range.mat
    }
  ### End internal functions
  
  if(!class(res.pca)[1]=="rda") stop("The input file is not a vegan output object of class 'rda'", call.=FALSE)
  if(scaling!=1 & scaling!=2) stop("Function only available for scaling = 1 or 2", call.=FALSE)
  
  k <- length(res.pca$CA$eig)         # n. of PCA eigenvalues
  n.sp <- length(res.pca$colsum)      # n. of species
  ahead <- 0.05   # Length of arrow heads
  aangle <- 30    # Angle of arrow heads
  # 'vec' will contain the selection of species to be drawn
  if(is.null(select.spe)){ vec <- 1:n.sp } else { vec <- select.spe }
  
  # Scaling 1: the species scores have norms of 1
  # Scaling 1: the site scores are scaled to variances = can.eigenvalues
  # Scaling 2: the species scores have norms of sqrt(can.eigenvalues)
  # Scaling 2: the site scores are scaled to variances of 1
  # --------------------------------------------------------------------
  
  ### This version reconstructs and uses the original RDA output of L&L 2012, Section 11.1.3
  
  Tot.var = res.pca$tot.chi        # Total variance in response data Y
  eig.val = res.pca$CA$eig         # Eigenvalues of Y-hat
  Lambda = diag(eig.val)           # Diagonal matrix of eigenvalues
  eig.val.rel = eig.val/Tot.var    # Relative eigenvalues of Y-hat
  Diag = diag(sqrt(eig.val.rel))   # Diagonal matrix of sqrt(relative eigenvalues)
  U.sc1 = res.pca$CA$v             # Species scores, scaling=1
  U.sc2 = U.sc1 %*% Lambda^(0.5)   # Species scores, scaling=2
  n = nrow(res.pca$CA$u)           # Number of observations
  Z.sc2 = res.pca$CA$u*sqrt(n-1)   # Site scores, scaling=2
  Z.sc1 = Z.sc2 %*% Lambda^(0.5)   # Site scores, scaling=1
  #
  if(is.null(select.spe)){ vec <- 1:n.sp } else { vec <- select.spe }
  #
  if(scaling==1) {
    sit.sc <- Z.sc1
    spe.sc <- U.sc1[vec,]
  } else {          # For scaling=2
    sit.sc <- Z.sc2
    spe.sc <- U.sc2[vec,]
  }
  if(is.null(rownames(sit.sc))) rownames(sit.sc) <- paste("Site",1:n,sep="")
  if(is.null(rownames(spe.sc))) rownames(spe.sc) <- paste("Sp",1:n.sp,sep="")
  #
  fact.spe <- 1
  if(optimum) {
    fact.spe <- stretch(sit.sc[,1:k], spe.sc[,1:k], ax1, ax2, n, silent=silent)
  }
  if(!silent) cat("fac.spe =",fact.spe,"\n\n")
  spe.sc <- spe.sc*fact.spe*mult.spe
  #
  lim <- larger.plot(sit.sc[,1:k], spe.sc[,1:k], percent=mar.percent, move.origin=move.origin, ax1=ax1, ax2=ax2)
  if(!silent) print(lim)
  
  ### Drawing the biplot begins ###
  ###
  # Draw the main plot
  mat <- rbind(sit.sc[,1:k], spe.sc[,1:k])
  plot(mat[,c(ax1,ax2)], type="n", main=paste("Biplot PCA, escalamiento", scaling), xlim=c(lim[1,ax1], lim[2,ax1]), ylim=c(lim[1,ax2], lim[2,ax2]), 
       xlab=paste("PCA ",ax1), ylab=paste("PCA ",ax2), asp=1)
  abline(h=0, v=0, col="grey60")
  
  # Draw the site scores
  if(plot.sites) {
    points(sit.sc[,ax1], sit.sc[,ax2], pch=20)
    if(label.sites)
      text(sit.sc[,ax1], sit.sc[,ax2], labels = rownames(sit.sc), col="black", pos=pos.sites, cex=cex.char1)
  } else {
    if(label.sites)
      text(sit.sc[,ax1], sit.sc[,ax2], labels = rownames(sit.sc), col="black", pos=NULL, cex=cex.char1)
  }
  
  # Draw the species scores
  if(plot.spe) {
    arrows(0, 0, spe.sc[,ax1], spe.sc[,ax2], length=ahead, angle=aangle, col="red")
    if(label.spe)
      text(spe.sc[,ax1], spe.sc[,ax2], labels = rownames(spe.sc), col="red", pos=pos.spe, cex=cex.char1)
  } else {
    if(label.spe)
      text(spe.sc[,ax1], spe.sc[,ax2], labels = rownames(spe.sc), col="red", pos=NULL, cex=cex.char1)
  }
  
  # If scaling=1 draw circle of equilibrium contribution
  #  if(scaling==1 | (scaling==2 & circle2)){
  if(scaling==1){
    pcacircle(res.pca, mult.spe=mult.spe, fact.spe=fact.spe, silent=silent)
  }
}

"pcacircle" <- function (pca, mult.spe, fact.spe, silent=silent) 
{
  # Draws a circle of equilibrium contribution on a PCA plot 
  # generated from a vegan analysis.
  
  eigenv <- pca$CA$eig
  p <- length(eigenv)
  n <- nrow(pca$CA$u)
  tot <- sum(eigenv)
  radius <- (2/p)^0.5 * mult.spe * fact.spe
  symbols(0, 0, circles=radius, inches=FALSE, add=TRUE, fg=2)
  if(!silent) cat("\nSpecies arrows and the radius of the equilibrium circle are stretched by a factor of", mult.spe*fact.spe)
  if(!silent) cat("\nThe radius of the equilibrium circle is thus", (2/p)^0.5, "*", mult.spe, "*", fact.spe, "=", radius,"\n")
}

dind <- function(mc, nomsitios=NA){
  require(vegan)
  require(BiodiversityR)
  require(fossil)
  zz1 <- cbind(
    ID=rownames(mc),
    as.data.frame(
      sapply(
        c("richness",
          "abundance",
          "Shannon",
          "Simpson",
          "inverseSimpson",
          "Logalpha",
          "Berger",
          "Jevenness",
          "Eevenness" ),
        function(x)
          if(
            package_version(packageVersion("BiodiversityR")) >
            package_version("2.9-1"))
            diversityresult(mc,index=x,method='each site', digits = 4) else
              diversityresult(mc,index=x,method='s', digits = 4)
      )
    )
  )
  chao <- data.frame(ID=rownames(mc), `estimador Chao1`=sapply(rownames(mc), function(x) chao1(mc[x,], taxa.row = F)))
  jack <- data.frame(ID=rownames(mc), `estimador Jackknife`=sapply(rownames(mc), function(x) jack1(mc[x,], taxa.row = F)))
  zz2 <- Reduce(function(x, y) merge(x, y, all=TRUE), list(zz1, chao, jack))
  # zz2 <- base::merge(zz, chao)
  colnames(zz2)[2:ncol(zz2)] <- c('riqueza','abundancia','Shannon','Gini-Simpson','inverso de Simpson','Fisher-alpha o Logalpha','Berger-Parker','equidad de Pielou o J-evenness','equidad de Buzas-Gibson o E-evenness', 'estimador Chao1', 'estimador Jackknife')
  if(!is.na(nomsitios)) colnames(zz2)[1] <- nomsitios
  return(zz2)
}

alpha_div <- function(mc = mi_fam) {
  # Tomado de: Borcard et al., 2018
  library(vegan)
  N0 <- specnumber(mc) # Riqueza de especies
  H <- diversity(mc) # Entropía de Shanon, base e
  Hb2 <- diversity(mc, base = 2) # Entropía de Shanon, base 2
  N1 <- exp(H) # Número de diversidad Hill 1, base e. Num. de sp. abundantes
  N1b2 <- 2^Hb2 # Número de diversidad de Hill 1, base 2.
  N2 <- diversity(mc, "inv") # Inverso de Simpson, diversidad. Num. de sp. dominantes
  J <- H / log(N0) # Equidad de Pielou
  E10 <- N1 / N0 # Equidad de Shannon (ratio de Hill)
  E20 <- N2 / N0 # Equidad de Simpson (ratio de Hill)
  div <- data.frame(N0, H, Hb2, N1, N1b2, N2, J, E10, E20)
  return(div)
}

estimacion_riqueza_chao <- function(mc, tamano_rarefaccion, n_raras = 10) {
  library(RColorBrewer)
  library(SpadeR)
  library(iNEXT)
  mc_lista <- sapply(
    rownames(mc), function(x) as.numeric(mc[x,]), simplify = F)
  asin_resumen_estimadores <- if(is.list(mc))
    sapply(
      mc_lista,
      function(x) 
        SpadeR::ChaoSpecies(x, datatype = 'abundance', k = n_raras, conf = 0.95),
      simplify = F)
  else
    SpadeR::ChaoSpecies(mc, datatype = 'abundance', k = n_raras, conf = 0.95)
   nasin_raref <- iNEXT::iNEXT(
    if(is.list(mc)) mc_lista else mc,
    q=0,
    knots = 400,
    datatype ="abundance",
    endpoint = tamano_rarefaccion)
  acumulacion_especies <- iNEXT::ggiNEXT(nasin_raref, type=1) +
    theme_bw() +
    theme(
      text = element_text(size = 20),
      panel.background = element_rect(fill = 'white', colour = 'black'),
      panel.grid.major = element_line(colour = "grey", linetype = "dashed", size = 0.25)
    ) +
    ylab('Riqueza de especies') +
    xlab('Número de individuos') +
    scale_y_continuous(breaks = seq(0,80, length.out = 9)) +
    scale_color_manual(values = brewer.pal(8, 'Set2')) +
    scale_fill_manual(values = brewer.pal(8, 'Set2'))
  return(list(
    asintoticos_estimacion = asin_resumen_estimadores,
    no_asintoticos_rarefaccion_extrapolacion = nasin_raref,
    no_asintoticos_rarefaccion_extrapolacion_grafico = acumulacion_especies))
}
