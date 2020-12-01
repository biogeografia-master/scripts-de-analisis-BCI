#' ---
#' title: "Análisis de diversidad. <br> Parte 1: Diversidad alpha"
#' author: "JR"
#' date: "29 de noviembre, 2020"
#' output: github_document
#' ---

knitr::opts_chunk$set(fig.width=12, fig.height=8)

#' ## Preámbulo
#' 
#' ### Cargar paquetes
#' 
library(vegan)
library(adespatial)
library(plyr)
library(RColorBrewer)
library(tidyverse)
library(sf)
library(SpadeR)
library(iNEXT)
source('biodata/funciones.R')
#' 
#' ### Cargar datos
#' 
load('biodata/Apocynaceae-Meliaceae-Sapotaceae.Rdata')
load('biodata/matriz_ambiental.Rdata')
mi_fam <- mc_apcyn_melic_saptc
bci_env_grid %>% tibble
grupos_upgma_k2 <- readRDS('grupos_upgma_k2.RDS')
table(grupos_upgma_k2)
grupos_ward_k3 <- readRDS('grupos_ward_k3.RDS')
table(grupos_ward_k3)
#' 
#' ## Diversidad alpha
#' 
#' **La tentación de medir la diversidad mediante un único número, ha atrapado a muchos investigadores en el pasado**. En distintos momentos, sobre todo durante el siglo XX, se desarrollaron varios índices de diversidad, que hoy vemos en la bibliografía, a veces usados de manera indiscriminada o sin sentido ecológico.
#' 
#' **La riqueza de especies, una cifra de uso común en estudios ecológicos, es un elemento muy simple** de la medición biodiversidad. La diversidad va más allá de la riqueza de especies, puesto que tiene muchas dimensiones y puede medirse de múltiples maneras. Además, **hay otros tipos de diversidad, desde el genoma hasta el paisaje**, y en estudios de ecología de comunidades, está atrayendo cada vez más atención el análisis de la **diversidad funcional y la diversidad filogenética** (Borcard et al., 2018; Magurran, 2004).
#' 
#' El término **"diversidad biológica"** se asume que fue acuñado a principios de los 80, y **se atribuye a Lovejoy**. Magurran lo atribuye a Gerbilskii y Petrunkevitch (1955), pero con una acepción diferente a la actual. Whiteside y Harmsworth (1967) lo usan también tempranamente. Norse y otros (1986), la dividen por primera vez en genética, de especies y ecológica.
#' 
#' **Biodiversidad** es de factura más reciente, y es una contracción de "biológica-diversidad" (del inglés *biological diversity*). Se atribuye a **Walter G. Rosen** (1986), mientras planeaba un evento científico sobre el tema. Las actas del evento se publicaron como un **libro titulado "Biodiversidad", editado por Edward O. Wilson**.
#' 
#' La definición de biodiversidad del PNUMA es: "**variabilidad entre organismos vivos** de todos los medios, incluyendo terrestres, marinos y otros sistemas acuáticos, y los complejos ecológicos de los cuales forman parte. Esto incluye la diversidad intraespecífica, interespecífica y de ecosistemas". Se trata de una definición de sentido amplio. Harper y Hawsworth sugieren esos mismos tres niveles, con los adjetivos **"genética", "de organismos" y "ecológica"**.
#' 
#' **Hubbell** (2001) ofrece una definición más adaptada a la práctica actual y mucho más restringida: "biodiversidad es sinónimo de **riqueza de especies y de abundancia relativa de especies** en el espacio y en el tiempo". Magurran utiliza **"diversidad biológica" y "biodiversidad" como sinónimos**, y la define como "la **variedad y la abundancia** de especies en una unidad de estudio", siguiendo lo planteado por Hubbel.
#' 
#' En estas dos últimas acepciones (Hubbell y Magurran), la diversidad biológica puede dividirse en dos componentes: **riqueza de especies y equidad**. Las mediciones de la biodiversidad, de las cuales hay un gran número, dan un peso relativo diferenciado a dichos componentes. 
#' 
#' ### La diversidad de especies como un único número
#' 
#' Usaré la notación *q* para designar el número de especies o riqueza de especies. Cualquier unidad de muestreo contiene un número determinado de individuos que pertenece a un cierto número de especies y, dado el hecho de que algunas especies son más raras que otras, es decir, son menos detectables, **el número total de especies de una unidad de muestreo o de un conjunto de unidades de muestreo, se incrementa al aumentar el área/volumen o el número de individuos muestreados**. Por lo tanto, la comparación de la riqueza de especies entre dos unidades de muestreo, la cual es un estimado del número de especies real, estará sesgada (Borcard et al., 2018).
#' 
#' #### Riqueza de especies y rarefacción
#' 
#' Magurran (2004) distingue entre **densidad de especies**, que equivale al número de especies por unidad de área de colecta, y **riqueza numérica de especies**, que es el número de especies por número de individuos o por unidad de biomasa.
#' 
#' Para asegurar la comparabilidad entre sitios, se han propuesto distintos métodos. Uno es la rarefacción, propuesta originalmente por Sanders (1968) y estandarizada por Hurlbert (1971), que estima el número de especies en unidades de muestreo conteniendo el mismo número de individuos, usando datos no transformados; se basa, por lo tanto, en el concepto de riqueza numérica de especies. Es decir, se determina *q'* por unidad estándar de muestreo *n'* de un universo que contiene *q* especies, *n* individuos y *n<sub>i</sub>* individuos pertenecientes a *i* especies.
#' 
#' ![](rarefaccion.jpg)
#' 
#' #### Componentes de la diversidad de especies basada en abundancia: riqueza y equidad
#' 
#' Si asumimos que un sitio de muestreo es una variable cualitativa, y cada especies un "estado". Bajo esta lógica, la dispersión de esta variable se calcula usando las frecuencias relativas *p<sub>i</sub>* de los *q*-estados usando la conocida entropía de Shannon (1948):
#' 
#' ![](shannon.jpg)
#' 
#' Desde el punto de vista ecológico, la entropía de Shannon tiene dos propiedades importantes: 1) Crece al aumentar la riqueza de especies *q*; 2) Crece con la uniformidad (=equidad o equitabilidad, es decir, qué tan bien repartida se encuentra la abundancia entre las especies). Para una *q* dada, la entropía de Shannon asume su valor máximo cuando todas las especies están igualmente representadas, y es equivalente al logaritmo de la riqueza:
#' 
#' ![](shannon_max.jpg)
#' 
foo1 <- c(25, 16, 9, 4, 1)
diversity(foo1)
foo2 <- c(11, 11, 11, 11, 11)
diversity(foo2)
#' 
#' Por otra parte, la **equidad de Pielou** (1966) es la razón entre la entropía de Shannon y su valor máximo; también se le conoce como **equidad de Shannon**. La equidad de Pielou es estrechamente dependiente de la riqueza, pero es un índice muy usado en trabajos ecológicos:
#'  
#' ![](pielou.jpg)
#'
#' La equidad (en general, no sólo la medida por Pielou), se relaciona con la forma de los modelos de abundancia de especies, que son funciones ajustadas a los gráficos a los gráficos rango-abundancia (horizontal=especies por rango de abundancia, vertical=logaritmo de las abundancias). Los principales modelos son, ordenados de menor a mayor equidad representada, geométrico, log, lognormal y de la vara quebrada. Los modelos de abundancia de especies se pueden consultar mediante la función `radfit` de `{vegan}`. La mayoría de estos modelos son realmente modelos lineales generalizados.
#' 
#' Otra medida común es el índice de concentración de Simpson (1949), λ, que equivale a la probabilidad de que dos individuos elegidos al azar pertenezcan a la misma especie.
#' 
#' ![](simpson.jpg)
#' 
#' Este valor aumenta con la dominancia (de ahí su nombre "índice de concentración"), por lo que realmente no mide diversidad, sino más bien dominancia. Para transformarlo en un índice de diversidad, se utiliza D=1-λ, que es el índice de Gini-Simpson, o D=1/λ, que es el inverso de Simpson (**menos sensible a cambios de la abundancia en las especies muy comunes**).
#' 
#' La riqueza de especies, la entropía de Shannon y la **diversidad** de Simpson son realmente casos especiales de la entropía generalizada de Renyi (1961):
#' 
#' ![](renyi.jpg)
#' 
#' donde *a* es el orden de la medida de entropía (*a=0,1,2...*), la cual cuantifica la importancia de la abundancia de especies y, por lo tanto, la equidad. Hill (1973) propuso usar los correspondientes números de diversidad:
#' 
#' ![](hill.jpg)
#' 
#' Las tres primeras entropias de Renyi (*H<sub>a</sub>*), donde *a=0, 1 y 2*, y los correspondientes números de diversidad de Hill, (*N<sub>a</sub>*), son realmente índices que ya conocemos: *H<sub>0</sub>=H<sub>max</sub>=log(q)*, *H<sub>1</sub>=H=entropia de Shannon*, *H<sub>2</sub>=-log(λ)*. Por otra parte, los tres primeros números de diversidad de Hill, *N<sub>0</sub>=q*, simplemente la riqueza de especies, *N<sub>1</sub>=exp(H)*, número de especies abundantes, y *N<sub>1</sub>=1/λ*, inverso de Simpson. De lo anterior se deriva que, **a medida que se incrementa *a*, se le da mayor importancia a la o las especies más abundantes**.
#' 
#' ![](tres_entro_renyi_hill_div_num.jpg) <br> 
#' > Según Borcard et al., 2018.
#' 
#' Bajo esta notación, el índice de equidad de Pielou (o equidad de Shannon) equivale a *J=H<sub>1</sub>/H<sub>0</sub>*, que es a fin de cuentas una ratio. Hill propuso también otras ratios: *E<sub>1</sub>=N<sub>1</sub>/N<sub>0</sub>* a la cual el propio Hill denominó como su versión de la equidad de Shannon y *E<sub>2</sub>=N<sub>2</sub>/N<sub>0</sub>*. Por lo tanto, Hill no sólo propuso números de diversidad, sino también ratios.
#' 
#' Los números de diversidad y las ratios de Hill son menos sensibles a las matrices de comunidad con fuerte dominancia, y producen los denominados "números equivalentes". Se pueden interpretar como "el número de elementos igualmente probables (individuos, especies, etc.) necesarios para producir el valor observado del índice de diversidad" (Ellison, 2010, modificado por Jost, 2007). Además, **los números de diversidad de Hill son preferibles para la interpretación a través de modelos lineales, porque tienen mayor probabilidad de estar relacionados linealmente con variables ambientales**.
#' 
#' **Estas afirmaciones tienen implicaciones muy importantes desde el punto de vista ecológico, puesto que tus datos podrían mostrar tendencias antes los números de Hill y no necesariamente con la entropia de Shannon o el clásico índice de Simpson**.
#' 
#' **Índices**
#' 
(indices <- alpha_div(mi_fam))
pairs(indices,
      lower.panel = panel.smooth,
      upper.panel = panel.cor,
      diag.panel = panel.hist,
      main = "Pearson Correlation Matrix")
indices_env <- bind_cols(
  indices,
  bci_env_grid %>%
    select_if(is.numeric) %>%
    st_drop_geometry %>%
    select(-id) %>% 
    select(-matches('^geom.*pct$')))
indices_env %>% tibble
ezCorM(indices_env, r_size_lims = c(3,5), label_size = 4)
#' 
#' **Modelos de abundancia de especies**
#' 
mi_fam_mae <- radfit(mi_fam)
plot(mi_fam_mae)
#' 
#' **Rarefacción**
#' 
#' Riqueza por sitio
#' 
riqueza <- specnumber(mi_fam)
riqueza %>% sort
#' 
#' Sitios con riqueza máxima y mínima
#' 
riqueza[riqueza == min(riqueza)]
riqueza[riqueza == max(riqueza)]
range(riqueza)
#' 
#' Abundancia por sitio
#' 
abundancia <- rowSums(mi_fam)
abundancia %>% sort
#' 
#' Sitios con abundancia máxima y mínima
#' 
abundancia[abundancia == min(abundancia)]
abundancia[abundancia == max(abundancia)]
(rango_abun <- range(abundancia))
#'
#' Abundancia en el sitio más pobre
#' 
abundancia[riqueza == min(riqueza)]
#' 
#' Abundancia en el sitio más rico
#' 
abundancia[riqueza == max(riqueza)]
#' 
#' Riqueza en el sitio con menor abundancia
#' 
riqueza[abundancia == min(abundancia)]
#' 
#' Riqueza en el sitio con mayor abundancia
#' 
riqueza[abundancia == max(abundancia)]
#' 
#' Rarefacción a la abundancia más pequeña encontrada
#' 
riqueza_menor_abun <- rarefy(mi_fam, sample = rango_abun[1])
# Compare ranking of observed and rarefied cores
sort(riqueza)
sort(round(riqueza_menor_abun))
rarecurve(
  mi_fam,
  step = 1,
  sample = rango_abun[1],
  xlab = "Número de individuos (tamaño de muestra)",
  ylab = "Especies",
  label = TRUE,
  col = "blue"
)
#' 
#' ### Riqueza de especies, estimación y comparación, "completitud de muestra" (existe en el diccionario) (Chao y Chiu, 2016)
#' 
specpool(mi_fam)
specpool(mi_fam)[[1]]/specpool(mi_fam)*100
#' 
#' Lista comprensiva de métodos (incluyendo recientes):
#' 
#' - **Enfoques asintóticos. Estiman la riqueza de especies**:
#'     - Paramétricos:
#'         - Modelo homogéneo (estándar y MLE)
#'     - No paramétricos:
#'         - Chao1 y Chao1-bc
#'         - iChao1
#'         - Basados en "cobertura" o "completitud de muestra". ACE para datos de abundancia
#'         - Estimadores Jackknife (de primer y segundo órdenes)
#' - **Enfoques no asintóticos. Se utilizan para hacer rarefacción y extrapolación**:
#'     - Basados en tamaño de la muestra
#'     - Basados en "cobertura" o "completitud de muestra"
#' 
#' Matriz de comunidad combinada (todos los sitios forman uno)
#' 
mi_fam_combinada <- colSums(mi_fam)
mi_fam_combinada %>% sort
mi_fam_combinada_chao <- estimacion_riqueza_chao(
  mc = mi_fam_combinada,
  tamano_rarefaccion = 40000)
mi_fam_combinada_chao$asintoticos_estimacion
mi_fam_combinada_chao$no_asintoticos_rarefaccion_extrapolacion
mi_fam_combinada_chao$no_asintoticos_rarefaccion_extrapolacion_grafico
#'
#' Matriz de comunidad agrupada según Ward (tres grupos)
#' 
mi_fam_k3 <- mi_fam %>%
  mutate(g=grupos_ward_k3) %>%
  group_by(g) %>%
  summarise_all(sum) %>%
  select(-g) %>% 
  data.frame
mi_fam_k3 %>% rowSums %>% sort
mi_fam_k3_chao <- estimacion_riqueza_chao(
  mc = mi_fam_k3,
  tamano_rarefaccion = 20000)
mi_fam_k3_chao$asintoticos_estimacion
mi_fam_k3_chao$no_asintoticos_rarefaccion_extrapolacion
mi_fam_k3_chao$no_asintoticos_rarefaccion_extrapolacion_grafico
