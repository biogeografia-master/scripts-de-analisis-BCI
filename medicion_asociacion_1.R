#' ---
#' title: "Medición de asociación"
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
load('biodata/matriz_ambiental.Rdata')
load('biodata/Apocynaceae-Meliaceae-Sapotaceae.Rdata')

#' ## Principales categorías de medidas de asociación: modos Q y R
#' 
#' En modo Q mides asociación entre pares de objetos, como por ejemplo, dos sitios de muestreo. Por lo tanto, en este caso, mides la asociación por medio de **la disimilaridad y la similaridad** entre pares de objetos, usando métricas como la **distancia euclídea o la similaridad de Jaccard**.
#' 
#' En cambio, en modo R, mides asociación entre pares de descriptores, como por ejemplo, dos variables, o dos especies (como variables). En este caso mides la asociación por medio de **la dependencia entre variables**, usando por ejemplo la **covarianza o el índice de correlación**.
#' 
#' ## Modo Q: matrices de disimilaridad entre objetos
#' 
#' ### Modo Q para datos cuantitativos de especies (abundancia)
#' 
#' __Las matrices de disimilaridad son matrices de distancia__. La disimilaridad (D) también la puedes obtener a partir de la similaridad (S), aplicando la fórmula D = 1 - S, y viceversa, S = 1 - D.
#' 
#' #### La paradoja de Orlóci
(mc_orloci <- tibble(
  sp1 = c(0, 1, 0),
  sp2 = c(1, 0, 4),
  sp3 = c(1, 0, 8),
  sitio = paste0('sit', 1:3)) %>% 
    column_to_rownames('sitio'))
#' Las matrices de distancia normalmente se muestran de la siguiente manera:
(dist.ldc(mc_orloci, "euclidean", silent = T))

#' Sin embargo, para facilitar la lectura de las distancias, en esta explicación ordeneré las matrices de distancia en columnas, usando la función de ayuda `organizar_matriz_distancia`. La primera que generaré es la de distancias euclideas a partir de datos brutos:
(d_euc <- dist.ldc(mc_orloci, "euclidean", silent = T) %>%
  organizar_matriz_distancia(func_dist = 'Euclidean'))

#' Siendo los sitios 1 y 2 tan diferentes en cuanto a las especies que los componen (no comparten especies de hecho), ¿por qué están tan próximos? Asimismo, los sitios 1 y 3 comparten especies, entonces, ¿por qué están tan distantes (=disímiles)? La explicación se atribuye a los ceros y a los valores de abundancia extremos. El hecho de que un par de sitios registren varias ausencias (o pseudo-ausencias), hace que "aparezcan" muy próximos (similares) en el espacio euclídeo. Esta paradoja sugiere que es necesario evitar la distancia euclídea a partir de datos brutos como métrica para comparar sitios.
#' 
#' Existen distintas maneras de solucionar este problema, normalmente recurriendo a métodos de estandarización de los datos brutos de la matriz de comunidad, y luego calculando distancia euclídea. Es decir, se obtiene una matriz transformada y a partir de ella se obtienen distancias. Los más comunes son el método de cuerdas (*chord*), *ji*-cuadrado y *Hellinger*. La función `dist.ldc` del paquete `adespatial` provee una solución en un único paso para cada caso.
#' 
#' - *Chord*:
d_cho <- dist.ldc(mc_orloci, "chord", silent = T) %>%
  organizar_matriz_distancia(func_dist = 'Chord')
#' - *Ji*-cuadrado:
d_chi <- dist.ldc(mc_orloci, "chisquare", silent = T) %>%
  organizar_matriz_distancia(func_dist = 'chi-square distance')
#' - *Hellinger* (valores primero divididos por abundancia total > sqrt)
d_hel <- dist.ldc(mc_orloci, "hellinger", silent = T) %>%
  organizar_matriz_distancia(func_dist = 'Hellinger')
#' - Uniendo y comparando
(d_todas <- bind_rows(d_euc, d_cho, d_chi, d_hel))
mc_orloci

#' Ahora observa cómo se calcula *chord* en dos pasos:
#' 
#' - Primero construye matriz de comunidad estandarizada, obtenida a partir del cuadrado de cada valor y luego dividido (estandarizando) por la suma de los cuadrados de toda la fila o sitio (vector unitario). El resultado final, es decir, la matriz con valores transformados, se obtiene a partir de la raíz cuadrada de cada valor.
(mc_orloci_norm <- sqrt(mc_orloci^2/rowSums(mc_orloci^2)))
#' La matriz de comunidad se dice que está "normalizada". Lo anterior se puede hacer más fácilmente con `decostand` de `vegan`:
(mc_orloci_norm <- decostand(mc_orloci, "normalize"))
#' - Luego, se obtiene la distancia euclidea:
(d_cho_2_pasos <- dist(mc_orloci_norm, method = 'euclidean')) %>% 
  organizar_matriz_distancia(func_dist = 'Chord en dos pasos')
d_cho #Compara

#' Aplicado a BCI y mi familia (en forma de matriz de distncia), utilizando la transformación *Hellinger*:
mi_fam_d_hel <- dist.ldc(mc_apcyn_melic_saptc, "hellinger", silent = T)

#' Para interpretar esta matriz, es necesario representarla esta matriz gráficamente:
coldiss(mi_fam_d_hel, diag = T)
#' Guardo el gráfico. Puedes usar el botón `Export` de la pestaña `Plots`
#' 
#' La calidad de gráficos exigida en revistas, suele requerir usar funciones específicas para guardar gráficos. Por ejemplo
png(
  filename = 'matriz_disimilaridad_hellinger.png',
  width = 2400, height = 1200, pointsize = 32
)
coldiss(mi_fam_d_hel, diag = T)
dev.off()
#' MUY IMPORTANTE. La última función, `dev.off()`, es necesaria para cerrar el dispositivo. Si no la ejecutas, no se generarán gráficos en el dispositivo estándar (pestaña `Plots`) 
#' 
#' ## Modo Q para datos binarios
#' 
mi_fam_d_pa <- vegdist(mc_apcyn_melic_saptc, method = 'jac', binary = T)
#' `binary=T` realiza primero `decostand(mc_apcyn_melic_saptc, method = 'pa')`
#' En esta matriz, valor pequeño significa que los sitios comparados son muy parecidos.
coldiss(mi_fam_d_pa, diag = T)
#














#' ### Asociación entre especies
#' 
mc_apcyn_melic_saptc_t <- t(mc_apcyn_melic_saptc)
mc_apcyn_melic_saptc_t_chi <- decostand(mc_apcyn_melic_saptc_t, "chi.square")
mc_apcyn_melic_saptc_t_chi_euc <- dist(mc_apcyn_melic_saptc_t_chi)
coldiss(mc_apcyn_melic_saptc_t_chi_euc, diag = TRUE)

#' ### Riqueza capturada según estimadores
specpool(mc_apcyn_melic_saptc)
specpool(mc_apcyn_melic_saptc)[[1]]/specpool(mc_apcyn_melic_saptc)*100

#SCBD (species contribution to beta diversity) y LCBD (local contribution...)
mc_apcyn_melic_saptc_beta <- beta.div(mc_apcyn_melic_saptc, method = "hellinger", nperm = 9999)
mc_apcyn_melic_saptc_beta$SCBD[mc_apcyn_melic_saptc_beta$SCBD >= mean(mc_apcyn_melic_saptc_beta$SCBD)]
row.names(mc_apcyn_melic_saptc[which(mc_apcyn_melic_saptc_beta$p.LCBD <= 0.05),])
#