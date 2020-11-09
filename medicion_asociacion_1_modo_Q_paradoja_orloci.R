#' ---
#' title: "Medición de asociación. Introducción a los modos de análisis Q y R. Modo Q aplicado a la paradoja de Orlóci"
#' author: "JR"
#' date: "3 de noviembre, 2020"
#' output: github_document
#' ---
#'
knitr::opts_chunk$set(fig.width=8, fig.height=5)
#'
#'  ## Preámbulo
#'  
#' ### Cargar paquetes
library(vegan)
library(adespatial)
library(tidyverse)
library(gridExtra)
source('biodata/funciones.R')
#' 
#' ## Modos Q y R
#' 
#' En modo Q mides asociación entre pares de objetos, como por ejemplo, entre dos sitios de muestreo. En este modo, mides la asociación por medio de **la disimilaridad y la similaridad** entre pares de objetos, usando métricas como la **distancia euclídea o la similaridad de Jaccard**.
#' 
#' En modo R mides asociación entre pares de descriptores, como por ejemplo, entre sdos variables, o dos especies. En este caso mides la asociación por medio de **la dependencia entre variables**, usando por ejemplo la **covarianza o el índice de correlación**.
#' 
#' ## Modo Q: matrices de disimilaridad entre objetos
#' 
#' ### Modo Q para datos cuantitativos de especies (abundancia). La paradoja de Orlóci
#' 
#' La paradoja de Orlóci (1978) plantea que la distancia euclidea es más pequeña entre dos sitios que no comparten especies que entre sitios que sí las comparten.
#' 
#' Esta paradoja se explica por la presencia de "ceros" (especies ausentes) en la matriz de comunidad, los cuales contribuyen a aumentar la similaridad de forma ficticia. Esto se comprueba de manera directa en la matriz de disimilaridad (=distancia). __La matriz de disimilaridad es simplemente una matriz de distancia__. La disimilaridad (D) también la puedes obtener a partir de la similaridad (S), aplicando la fórmula D = 1 - S, y viceversa, S = 1 - D.
#'
#' Te muestro la paradoja con un ejemplo y, posteriormente, te explico cómo solucionar el problema de los ceros. Sea una matriz de comunidad `mc_orloci` de 2 especies y tres sitios...
#' 
(mc_orloci <- tibble(
  sp1 = c(1, 0, 4),
  sp2 = c(1, 0, 8),
  sitio = paste0('sit', 1:3)) %>% 
    column_to_rownames('sitio'))
#' 
#' ...donde ambas especies están ausentes en `sit2`, en `sit1` presentes con poca abundancia y en `sit3` con valores relativamente extremos.
#' 
#' En modo Q, calcularé la "distancia" o "disimilaridad" entre sitios según las especies que los caracterizan. En el lenguaje R, es común usar la función `dist` para calcular la matriz de distancias, pero en ecología se usa también `dist.ldc` del paquete `adespatial`. Las matrices de distancia normalmente se muestran de la siguiente manera (por defecto, sólo el triángulo inferior):
#' 
(dist.ldc(mc_orloci, "euclidean", silent = T))
#' 
#' Te muestro un gráfico de dispersión de los sitios según la abundancia de especies (los ejes representan la abundancia de cada especie). Puedes comprobar que la distancia entre `sit1` y `sit2` es pequeña, mientras que entre `sit1` y `sit3` es grande.
#' 
mc_orloci %>% rownames_to_column('id') %>% 
  ggplot() +
  aes(x = sp1, y = sp2, label = id) +
  geom_point(size = 3) +
  geom_text(vjust="inward",hjust="inward", size = 5, color = 'grey40') +
  coord_equal() +
  theme_bw() +
  theme(text = element_text(size = 16))
#'
#' Para facilitar la lectura de las distancias, en esta explicación ordeneré las matrices de distancia en columnas, usando la función de ayuda `organizar_matriz_distancia`. La primera que generaré es la de distancias euclideas a partir de datos brutos:
#' 
(d_euc <- dist.ldc(mc_orloci, "euclidean", silent = T) %>%
  organizar_matriz_distancia(func_dist = 'Euclidean'))
#' 
#' Siendo los sitios 1 y 2 tan diferentes en cuanto a las especies que los componen (de hecho, no comparten especies), ¿por qué están tan próximos? Asimismo, los sitios 1 y 3 comparten especies, entonces, ¿por qué están tan distantes (=disímiles)? La explicación se atribuye a los ceros y a los valores de abundancia extremos. El hecho de que un par de sitios registren varias ausencias (o pseudo-ausencias), hace que "aparezcan" muy próximos (similares) en el espacio euclídeo. Esta paradoja sugiere que es necesario evitar la distancia euclídea a partir de datos brutos como métrica para comparar sitios.
#' 
#' Existen distintas maneras de solucionar el problema planteado en la paradoja, normalmente recurriendo a métodos de estandarización de las abundancia brutas, y luego calculando distancia euclídea. Es decir, se obtiene una matriz de comunidad transformada y, a partir de ella, se obtienen distancias. Los métodos de transformación más comunes son el de cuerdas (*chord*), *ji*-cuadrado y *Hellinger*. La función `dist.ldc` del paquete `adespatial` provee una solución en un único paso para cada caso (transformar + calcular distancias).
#' 
#' - *Chord*:
#' 
d_cho <- dist.ldc(mc_orloci, "chord", silent = T) %>%
  organizar_matriz_distancia(func_dist = 'Chord')
#' 
#' - *Ji*-cuadrado:
#' 
d_chi <- dist.ldc(mc_orloci, "chisquare", silent = T) %>%
  organizar_matriz_distancia(func_dist = 'chi-square distance')
#' 
#' - *Hellinger* (valores primero divididos por abundancia total > sqrt)
#' 
d_hel <- dist.ldc(mc_orloci, "hellinger", silent = T) %>%
  organizar_matriz_distancia(func_dist = 'Hellinger')
#' 
#' - Uniendo y comparando
(d_todas <- bind_rows(d_euc, d_cho, d_chi, d_hel))
#' 
#' Verás que el par `sit1|sit3` tiene corta distancia, es decir, son muy parecidos (0.17 en *Hellinger*, e igualmente, corta distancia relativa en *ji*-cuadrado y en *chord*). Esto se debe a la transformación, basada en estandarización. Es lógico que este par sea muy similar, puesto que comparten las únicas dos especies de la comunidad. Por lo tanto, las tres transformaciones empleadas corrigen el problema planteado en la paradoja de Orlóci.
#' 
#' Nota igualmente que, tanto los pares `sit1`|`sit2` y `sit2|sit3` están distantes (distancia 1), lo cual también es deseable. Por tal razón, snos interesa explorar las matrices de comunidad original y la transformada. La función `dist.ldc` realiza dos pasos: primero genera una matriz de comunidad transformada y luego calcula la distancia. Para replicar el proceso que realiza `dist.ldc`, es necesario generar primero la matriz transformada y luego pasar directamente a calcular distancia. Repasa primero la matriz de comunidad original:
#' 
mc_orloci
#' 
#' A continuación, generaré la matriz transformada según el método *chord*. Esta matriz se calcula obteniendo primero el cuadrado de cada valor de abundancia y luego dividiéndolo (estandarizando) por la suma de los cuadrados de toda la fila o sitio (vector unitario). El resultado final, es decir, la matriz con valores transformados, se obtiene a partir de la raíz cuadrada de cada valor.
#' 
mc_orloci_norm <- sqrt(mc_orloci^2/rowSums(mc_orloci^2)) %>%
  replace(is.na(.), 0)
mc_orloci_norm
#' 
#' La matriz de comunidad se dice que está "normalizada". Lo anterior se puede hacer más fácilmente con `decostand` de `vegan` (resultado idéntico):
#' 
(mc_orloci_norm <- decostand(mc_orloci, "normalize"))
#' 
#' Al graficar los sitios sobre un espacio bidimensional, cada eje representando una especie, se obtienen resultados diferentes a los que se obtuvieron con la matriz original:
#' 
p1 <- mc_orloci %>%
  rownames_to_column('id') %>% 
  ggplot() +
  aes(x = sp1, y = sp2, label = id) +
  geom_point(size = 3) +
  geom_text(vjust="inward",hjust="inward", size = 5, color = 'grey40') +
  coord_equal() +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ggtitle('mc original')
p2 <- mc_orloci_norm %>%
  rownames_to_column('id') %>% 
  ggplot() +
  aes(x = sp1, y = sp2, label = id) +
  geom_point(size = 3) +
  geom_text(vjust="inward",hjust="inward", size = 5, color = 'grey40') +
  coord_equal() +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  ggtitle('mc transformada')
grid.arrange(p1, p2, nrow = 1)
#'   
#' - Por último, para completar el proceso realizado por `dist.ldc`, debes calcular la distancia euclidea a partir de esta matriz. Verás que obtienes el mismo resultado que con la función `dist.ldc`:
#' 
(d_cho_2_pasos <- dist(mc_orloci_norm, method = 'euclidean')) %>% 
  organizar_matriz_distancia(func_dist = 'Chord en dos pasos')
#' Compara la anterior con la generada por `dist.ldc`:
d_cho
#' 