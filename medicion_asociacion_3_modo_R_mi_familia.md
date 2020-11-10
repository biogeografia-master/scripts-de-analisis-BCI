Medición de asociación. Modo R aplicado a mi familia asignada
================
JR
3 de noviembre, 2020

``` r
knitr::opts_chunk$set(fig.width=12, fig.height=8)
```

## Preámbulo

### Cargar paquetes

``` r
library(vegan)
```

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.5-6

``` r
library(adespatial)
```

    ## Registered S3 methods overwritten by 'adegraphics':
    ##   method         from
    ##   biplot.dudi    ade4
    ##   kplot.foucart  ade4
    ##   kplot.mcoa     ade4
    ##   kplot.mfa      ade4
    ##   kplot.pta      ade4
    ##   kplot.sepan    ade4
    ##   kplot.statis   ade4
    ##   scatter.coa    ade4
    ##   scatter.dudi   ade4
    ##   scatter.nipals ade4
    ##   scatter.pco    ade4
    ##   score.acm      ade4
    ##   score.mix      ade4
    ##   score.pca      ade4
    ##   screeplot.dudi ade4

    ## Registered S3 method overwritten by 'spdep':
    ##   method   from
    ##   plot.mst ape

    ## Registered S3 methods overwritten by 'adespatial':
    ##   method             from       
    ##   plot.multispati    adegraphics
    ##   print.multispati   ade4       
    ##   summary.multispati ade4

``` r
library(broom)
library(tidyverse)
```

    ## ── Attaching packages ────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✓ ggplot2 3.3.2     ✓ purrr   0.3.4
    ## ✓ tibble  3.0.3     ✓ dplyr   0.8.3
    ## ✓ tidyr   1.0.0     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.4.0

    ## ── Conflicts ───────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(sf)
```

    ## Linking to GEOS 3.6.2, GDAL 2.2.3, PROJ 4.9.3

``` r
library(gclus)
```

    ## Loading required package: cluster

    ## Registered S3 method overwritten by 'gclus':
    ##   method         from 
    ##   reorder.hclust vegan

``` r
source('biodata/funciones.R')
```

### Cargar datos

``` r
load('biodata/matriz_ambiental.Rdata')
load('biodata/Apocynaceae-Meliaceae-Sapotaceae.Rdata')
```

## Modo R: matrices de dependencia entre variables (índice de correlación)

### Modo R para datos cuantitativos de especies (abundancia)

En este caso, las variables usaré los valores de abundancias de especies
como variables. Es decir, **compararé el grado de asociación entre
especies, NO entre sitios**.

Aunque se podría usar el índice de correlación como métrica de la
dependencia (tal como mostré en el script
`aed_5_correlaciones_variables_ambientales.R`), tanto los doble-ceros
(ausencias de una misma especie en dos lugares), como los *outliers* de
abundancia, contribuirían a aumentar de manera ficticia el índice de
correlación de Pearson.

Por tal razón, es recomendable aplicar la transformación *Chi* a la
matriz de comunidad transpuesta. Al utilizar una matriz transpuesto,
lograrás comparar especies, NO sitios (recuerda que en modo R comparamos
descriptores, no objetos). Explico el procedimiento a continuación, paso
a paso.

Primero, sustituyo el caracter de espacio por un <enter> en los nombres
de las especies (caracter \\n), para facilitar la lectura de los nombres
en la diagonal del mapa de calor. Luego transpongo la matriz.

``` r
mi_fam_t <- mc_apcyn_melic_saptc %>% 
  rename_all(gsub, pattern = ' ', replacement = '\n') %>% 
  t()
mi_fam_t %>% tibble
```

    ## # A tibble: 16 x 1
    ##    .[,"1"] [,"2"] [,"3"] [,"4"] [,"5"] [,"6"] [,"7"] [,"8"] [,"9"] [,"10"]
    ##      <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl>
    ##  1       3      2      2      3      2      3      0      4      5       1
    ##  2       0      0      0      0      1      0      0      0      0       0
    ##  3      21     11     19     38     21     18      9     14     25      17
    ##  4       2      1      1      2      0      0      3      0      0       1
    ##  5      24     20     15     19     24     12      8     10      5      19
    ##  6       4      2      2      1      1      1      0      1      1       0
    ##  7      37     33     20     21     21     15     17     24      5      17
    ##  8       2      0      3      3      3      3      2      4      1       1
    ##  9       0      0      0      0      1      0      0      0      0       0
    ## 10      12     23     29     18     16     19     31     17      7      25
    ## 11       0      0      1      4      1      2      1      3      1       2
    ## 12       0      0      0      0      0      0      0      0      0       0
    ## 13      37     37     51     48     41     43     36     37     26      42
    ## 14       6      1      0      0      2      0      1      0      0       0
    ## 15       2      6      5      8      7      2     14      7      9       4
    ## 16     159    177     67     72     47     75    131     94     95      85
    ## # … with 40 more variables: [,"11"] <dbl>, [,"12"] <dbl>, [,"13"] <dbl>,
    ## #   [,"14"] <dbl>, [,"15"] <dbl>, [,"16"] <dbl>, [,"17"] <dbl>,
    ## #   [,"18"] <dbl>, [,"19"] <dbl>, [,"20"] <dbl>, [,"21"] <dbl>,
    ## #   [,"22"] <dbl>, [,"23"] <dbl>, [,"24"] <dbl>, [,"25"] <dbl>,
    ## #   [,"26"] <dbl>, [,"27"] <dbl>, [,"28"] <dbl>, [,"29"] <dbl>,
    ## #   [,"30"] <dbl>, [,"31"] <dbl>, [,"32"] <dbl>, [,"33"] <dbl>,
    ## #   [,"34"] <dbl>, [,"35"] <dbl>, [,"36"] <dbl>, [,"37"] <dbl>,
    ## #   [,"38"] <dbl>, [,"39"] <dbl>, [,"40"] <dbl>, [,"41"] <dbl>,
    ## #   [,"42"] <dbl>, [,"43"] <dbl>, [,"44"] <dbl>, [,"45"] <dbl>,
    ## #   [,"46"] <dbl>, [,"47"] <dbl>, [,"48"] <dbl>, [,"49"] <dbl>,
    ## #   [,"50"] <dbl>

Segundo, transformo la matriz transpuesta usando estandarización *Chi*.

``` r
mi_fam_t_chi <- decostand(mi_fam_t, "chi.square")
mi_fam_t_chi %>% tibble
```

    ## # A tibble: 16 x 1
    ##    .[,"1"] [,"2"] [,"3"] [,"4"] [,"5"] [,"6"] [,"7"] [,"8"] [,"9"] [,"10"]
    ##      <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl>
    ##  1  0.0490 0.0324 0.0391 0.0559 0.0419 0.0620 0      0.0783 0.107   0.0196
    ##  2  0      0      0      0      0.825  0      0      0      0       0     
    ##  3  0.228  0.119  0.247  0.471  0.292  0.247  0.108  0.182  0.356   0.222 
    ##  4  0.0903 0.0449 0.0541 0.103  0      0      0.150  0      0       0.0543
    ##  5  0.256  0.212  0.192  0.231  0.328  0.162  0.0942 0.128  0.0698  0.243 
    ##  6  0.475  0.236  0.285  0.136  0.152  0.150  0      0.142  0.156   0     
    ##  7  0.151  0.134  0.0980 0.0980 0.110  0.0776 0.0768 0.118  0.0268  0.0835
    ##  8  0.151  0      0.272  0.259  0.291  0.287  0.167  0.363  0.0992  0.0910
    ##  9  0      0      0      0      3.30   0      0      0      0       0     
    ## 10  0.0855 0.163  0.248  0.146  0.146  0.171  0.244  0.145  0.0653  0.214 
    ## 11  0      0      0.154  0.588  0.165  0.326  0.142  0.463  0.169   0.309 
    ## 12  0      0      0      0      0      0      0      0      0       0     
    ## 13  0.165  0.164  0.273  0.244  0.234  0.243  0.177  0.198  0.152   0.225 
    ## 14  0.552  0.0913 0      0      0.236  0      0.102  0      0       0     
    ## 15  0.0327 0.0975 0.0981 0.149  0.147  0.0414 0.253  0.137  0.193   0.0786
    ## 16  0.113  0.125  0.0572 0.0586 0.0429 0.0676 0.103  0.0803 0.0887  0.0727
    ## # … with 40 more variables: [,"11"] <dbl>, [,"12"] <dbl>, [,"13"] <dbl>,
    ## #   [,"14"] <dbl>, [,"15"] <dbl>, [,"16"] <dbl>, [,"17"] <dbl>,
    ## #   [,"18"] <dbl>, [,"19"] <dbl>, [,"20"] <dbl>, [,"21"] <dbl>,
    ## #   [,"22"] <dbl>, [,"23"] <dbl>, [,"24"] <dbl>, [,"25"] <dbl>,
    ## #   [,"26"] <dbl>, [,"27"] <dbl>, [,"28"] <dbl>, [,"29"] <dbl>,
    ## #   [,"30"] <dbl>, [,"31"] <dbl>, [,"32"] <dbl>, [,"33"] <dbl>,
    ## #   [,"34"] <dbl>, [,"35"] <dbl>, [,"36"] <dbl>, [,"37"] <dbl>,
    ## #   [,"38"] <dbl>, [,"39"] <dbl>, [,"40"] <dbl>, [,"41"] <dbl>,
    ## #   [,"42"] <dbl>, [,"43"] <dbl>, [,"44"] <dbl>, [,"45"] <dbl>,
    ## #   [,"46"] <dbl>, [,"47"] <dbl>, [,"48"] <dbl>, [,"49"] <dbl>,
    ## #   [,"50"] <dbl>

Tercero, calculo la distancia euclídea.

``` r
mi_fam_t_chi_d <- dist(mi_fam_t_chi)
mi_fam_t_chi_d %>% tidy
```

    ## # A tibble: 120 x 3
    ##    item1                      item2                      distance
    ##    <fct>                      <fct>                         <dbl>
    ##  1 "Cedrela\nodorata"         "Aspidosperma\nspruceanum"     2.71
    ##  2 "Chrysophyllum\nargenteum" "Aspidosperma\nspruceanum"     1.21
    ##  3 "Chrysophyllum\ncainito"   "Aspidosperma\nspruceanum"     1.21
    ##  4 "Guarea\nbullata"          "Aspidosperma\nspruceanum"     1.10
    ##  5 "Guarea\ngrandifolia"      "Aspidosperma\nspruceanum"     1.40
    ##  6 "Guarea\nguidonia"         "Aspidosperma\nspruceanum"     1.11
    ##  7 "Lacmellea\npanamensis"    "Aspidosperma\nspruceanum"     1.42
    ##  8 "Pouteria\nfossicola"      "Aspidosperma\nspruceanum"     5.16
    ##  9 "Pouteria\nreticulata"     "Aspidosperma\nspruceanum"     1.08
    ## 10 "Pouteria\nstipitata"      "Aspidosperma\nspruceanum"     1.91
    ## # … with 110 more rows

Finalmente, creo el “mapa de
calor”.

``` r
coldiss(mi_fam_t_chi_d, diag = TRUE)
```

![](medicion_asociacion_3_modo_R_mi_familia_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

En el mapa de calor **ordenado** (el de la derecha), se identifica al
menos un patrón de dependencia entre las especies relacionadas en la
diagonal desde *Chrysophyllum cainito* hasta *Trichilia pallida*
(cuadros de color rosa centrales). También se observan las especies que
no parecen asociarse con otras, situadas en los extremos de la diagonal,
y relacionadas con otras por medio de valores pequeños de distancia
(cuadros azules), como *Rauvolfia littoralis* y *Pouteria fossicola* y
*Cedrela odorata*.

### Modo R para datos binarios (presencia/ausencia)

Arriba usé la distancia de Jaccard para evaluar asociación entre sitios.
Dicha métrica también se puede usar para evaluar la distancia entre
especies, usando como fuente la matriz de comunidad transpuesta
convertida a binaria (presencia/ausencia)

``` r
mi_fam_t_jac <- vegdist(mi_fam_t, "jaccard", binary = TRUE)
mi_fam_t_jac %>% tidy
```

    ## # A tibble: 120 x 3
    ##    item1                      item2                      distance
    ##    <fct>                      <fct>                         <dbl>
    ##  1 "Cedrela\nodorata"         "Aspidosperma\nspruceanum"   0.813 
    ##  2 "Chrysophyllum\nargenteum" "Aspidosperma\nspruceanum"   0.0400
    ##  3 "Chrysophyllum\ncainito"   "Aspidosperma\nspruceanum"   0.16  
    ##  4 "Guarea\nbullata"          "Aspidosperma\nspruceanum"   0.0400
    ##  5 "Guarea\ngrandifolia"      "Aspidosperma\nspruceanum"   0.271 
    ##  6 "Guarea\nguidonia"         "Aspidosperma\nspruceanum"   0.0400
    ##  7 "Lacmellea\npanamensis"    "Aspidosperma\nspruceanum"   0.220 
    ##  8 "Pouteria\nfossicola"      "Aspidosperma\nspruceanum"   0.938 
    ##  9 "Pouteria\nreticulata"     "Aspidosperma\nspruceanum"   0.0400
    ## 10 "Pouteria\nstipitata"      "Aspidosperma\nspruceanum"   0.420 
    ## # … with 110 more rows

``` r
coldiss(mi_fam_t_jac, diag = TRUE)
```

![](medicion_asociacion_3_modo_R_mi_familia_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

### Modo R para datos cuantitativos, NO de abundancia de especies (variables ambientales)

En modo R evalúas asociación entre descriptores, es decir, entre
variables. La métrica comúnmente usada es el índice de correlación de
Pearson. Sin embargo, si los datos no presentan distribución normal,
puedes emplear métricas más flexibles, como el índice *rho* de Spearman
o **tau** de Kendall.

En este ejemplo, mostraré la correlación entre variables de suelo y la
abundancia y riqueza globales y de mi familia asignada. Haré lo propio
con variables geomorfológicas. Ya tuviste ocasión de usar estos métodos
en el [*script* de análisis exploratorio de datos, sección
correlación](aed_5_correlaciones_variables_ambientales.md). En este
*script*, añadirás al análisis el índice *rho* de Spearman.

``` r
env_num <- bci_env_grid %>%
  dplyr::select_if(is.numeric) %>%
  dplyr::select(-id, -matches('^U.*')) %>% 
  st_drop_geometry %>% 
  mutate(
    riqueza_mifam = specnumber(mc_apcyn_melic_saptc),
    abundancia_mifam = rowSums(mc_apcyn_melic_saptc)) %>% 
  rename_all(gsub, pattern = '_pct$', replacement = '') %>% 
  rename_all(gsub, pattern = '_| ', replacement = '\n')
env_num %>% tibble
```

    ## # A tibble: 50 x 33
    ##    `heterogeneidad… `geomorf\nllanu… `geomorf\npico` `geomorf\ninter…
    ##               <dbl>            <dbl>           <dbl>            <dbl>
    ##  1           0.627             10.0             0                0.83
    ##  2           0.394             34.8             0                0.36
    ##  3           0                  0               0                0   
    ##  4           0                  0               0                0.16
    ##  5           0.461              2.58            0                0   
    ##  6           0.0768             0               0.17             3.01
    ##  7           0.381              0               0.53             2.87
    ##  8           0.211              0               0                0   
    ##  9           0                  0               0                0   
    ## 10           0                  1.03            0                0   
    ## # … with 40 more rows, and 29 more variables: `geomorf\nhombrera` <dbl>,
    ## #   `geomorf\nespolón/gajo` <dbl>, `geomorf\nvertiente` <dbl>,
    ## #   `geomorf\nvaguada` <dbl>, `geomorf\npiedemonte` <dbl>,
    ## #   `geomorf\nvalle` <dbl>, `geomorf\nsima` <dbl>, Al <dbl>, B <dbl>,
    ## #   Ca <dbl>, Cu <dbl>, Fe <dbl>, K <dbl>, Mg <dbl>, Mn <dbl>, P <dbl>,
    ## #   Zn <dbl>, N <dbl>, N.min. <dbl>, pH <dbl>, `elevacion\nmedia` <dbl>,
    ## #   `pendiente\nmedia` <dbl>, `orientacion\nmedia` <dbl>,
    ## #   `curvatura\nperfil\nmedia` <dbl>,
    ## #   `curvatura\ntangencial\nmedia` <dbl>, `abundancia\nglobal` <dbl>,
    ## #   `riqueza\nglobal` <int>, `riqueza\nmifam` <int>,
    ## #   `abundancia\nmifam` <dbl>

``` r
p_cor_suelo_ar <- env_num %>%
  dplyr::select(matches('^[A-T,Z]|abundancia|riqueza|^pH$', ignore.case = F)) %>%
  ezCorM(r_size_lims = c(4,8), label_size = 3, method = 'pearson')
```

    ## -------------------------------------------------------------------------

    ## You have loaded plyr after dplyr - this is likely to cause problems.
    ## If you need functions from both plyr and dplyr, please load plyr first, then dplyr:
    ## library(plyr); library(dplyr)

    ## -------------------------------------------------------------------------

    ## 
    ## Attaching package: 'plyr'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     arrange, count, desc, failwith, id, mutate, rename, summarise,
    ##     summarize

    ## The following object is masked from 'package:purrr':
    ## 
    ##     compact

    ## 
    ## Attaching package: 'reshape2'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     smiths

``` r
p_cor_suelo_ar
```

    ## `geom_smooth()` using formula 'y ~ x'

    ## `geom_smooth()` using formula 'y ~ x'

![](medicion_asociacion_3_modo_R_mi_familia_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
p_cor_suelo_ar_spearman <- env_num %>%
  dplyr::select(matches('^[A-T,Z]|abundancia|riqueza|^pH$', ignore.case = F)) %>%
  ezCorM(r_size_lims = c(4,8), label_size = 3, method = 'spearman')
p_cor_suelo_ar_spearman
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](medicion_asociacion_3_modo_R_mi_familia_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

``` r
p_cor_geomorf_ar <- env_num %>%
  dplyr::select(-matches('^[A-T,Z]|pH', ignore.case = F)) %>%
  ezCorM(r_size_lims = c(4,8), label_size = 3, method = 'pearson')
p_cor_geomorf_ar
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](medicion_asociacion_3_modo_R_mi_familia_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->

``` r
p_cor_geomorf_ar_spearman <- env_num %>%
  dplyr::select(-matches('^[A-T,Z]|pH', ignore.case = F)) %>%
  ezCorM(r_size_lims = c(4,8), label_size = 3, method = 'spearman')
p_cor_geomorf_ar_spearman
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](medicion_asociacion_3_modo_R_mi_familia_files/figure-gfm/unnamed-chunk-9-4.png)<!-- -->
