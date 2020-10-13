Análisis exploratorio de datos, parte 1
================
JR
13 de octubre, 2020

### Área de cargar paquetes

``` r
library(vegan)
```

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.5-6

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────── tidyverse 1.2.1 ──

    ## ✓ ggplot2 3.3.2     ✓ purrr   0.3.4
    ## ✓ tibble  3.0.3     ✓ dplyr   0.8.3
    ## ✓ tidyr   1.0.0     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.4.0

    ## ── Conflicts ────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(sf)
```

    ## Linking to GEOS 3.6.2, GDAL 2.2.3, PROJ 4.9.3

``` r
source('biodata/funciones.R')
```

### Área de cargar datos

Censo (el objeto se carga con prefijo “censo”) y matriz de comunidad
(prefijo “mc”)

``` r
load('biodata/Apocynaceae-Meliaceae-Sapotaceae.Rdata')
load('biodata/matriz_ambiental.Rdata') #Matriz ambiental, se carga como "bci_env_grid"
```

### Imprimir datos en pantalla (impresiones parciales con head)

``` r
head(censo_apcyn_melic_saptc)
```

    ##   treeID stemID    tag StemTag     sp quadrat    gx    gy MeasureID
    ## 1     35      1 000018         tri2tu    4921 991.6 431.8        12
    ## 2     41      4 000024    <NA> guargu    4921 991.4 420.2        18
    ## 3     47      1 000030         tri2tu    4919 987.9 391.5        27
    ## 4     75      1 000060         tri2tu    4915 993.8 307.5        47
    ## 5     90      1 000075         poutre    4912 987.0 250.4        58
    ## 6     99      1 000084         tab2ar    4912 982.1 240.7        62
    ##   CensusID dbh pom hom  ExactDate DFstatus  codes nostems  date status
    ## 1      171 398 1.3 1.3 2010-09-21    alive   <NA>       1 18526      A
    ## 2      171 287 1.3 1.3 2010-09-21    alive   <NA>       1 18526      A
    ## 3      171 328 1.3 1.3 2010-09-20    alive   <NA>       1 18525      A
    ## 4      171 378 3.2 3.3 2010-09-14    alive B,cylY       1 18519      A
    ## 5      171 545 1.3 1.3 2010-09-10    alive   <NA>       1 18515      A
    ## 6      171 338 3.2 3.2 2010-09-14    alive B,cylN       1 18519      A
    ##         agb                   Latin           Genus     Species
    ## 1 1.5844933   Trichilia tuberculata       Trichilia tuberculata
    ## 2 0.6108122         Guarea guidonia          Guarea    guidonia
    ## 3 0.9678957   Trichilia tuberculata       Trichilia tuberculata
    ## 4 1.3904099   Trichilia tuberculata       Trichilia tuberculata
    ## 5 3.8887012     Pouteria reticulata        Pouteria  reticulata
    ## 6 0.8604040 Tabernaemontana arborea Tabernaemontana     arborea
    ##        Family speciesID                 authority IDlevel
    ## 1   Meliaceae      1169 (Triana & Planch.) C. DC. species
    ## 2   Meliaceae       435              (L.) Sleumer species
    ## 3   Meliaceae      1169 (Triana & Planch.) C. DC. species
    ## 4   Meliaceae      1169 (Triana & Planch.) C. DC. species
    ## 5  Sapotaceae      1155              (Engl.) Eyma species
    ## 6 Apocynaceae       990         Rose in Donn. Sm. species
    ##                     syn subsp quad1ha
    ## 1        Trichilia cipo  <NA>      50
    ## 2                        <NA>      50
    ## 3        Trichilia cipo  <NA>      49
    ## 4        Trichilia cipo  <NA>      49
    ## 5 Pouteria unilocularis  <NA>      48
    ## 6                        <NA>      48

``` r
head(mc_apcyn_melic_saptc)
```

    ##   Aspidosperma spruceanum Cedrela odorata Chrysophyllum argenteum
    ## 1                       3               0                      21
    ## 2                       2               0                      11
    ## 3                       2               0                      19
    ## 4                       3               0                      38
    ## 5                       2               1                      21
    ## 6                       3               0                      18
    ##   Chrysophyllum cainito Guarea bullata Guarea grandifolia Guarea guidonia
    ## 1                     2             24                  4              37
    ## 2                     1             20                  2              33
    ## 3                     1             15                  2              20
    ## 4                     2             19                  1              21
    ## 5                     0             24                  1              21
    ## 6                     0             12                  1              15
    ##   Lacmellea panamensis Pouteria fossicola Pouteria reticulata
    ## 1                    2                  0                  12
    ## 2                    0                  0                  23
    ## 3                    3                  0                  29
    ## 4                    3                  0                  18
    ## 5                    3                  1                  16
    ## 6                    3                  0                  19
    ##   Pouteria stipitata Rauvolfia littoralis Tabernaemontana arborea
    ## 1                  0                    0                      37
    ## 2                  0                    0                      37
    ## 3                  1                    0                      51
    ## 4                  4                    0                      48
    ## 5                  1                    0                      41
    ## 6                  2                    0                      43
    ##   Thevetia ahouai Trichilia pallida Trichilia tuberculata
    ## 1               6                 2                   159
    ## 2               1                 6                   177
    ## 3               0                 5                    67
    ## 4               0                 8                    72
    ## 5               2                 7                    47
    ## 6               0                 2                    75

``` r
bci_env_grid # No necesita imprimirse parcialmente
```

    ## Simple feature collection with 50 features and 14 fields
    ## geometry type:  POLYGON
    ## dimension:      XY
    ## bbox:           xmin: 625704 ymin: 1011519 xmax: 626704 ymax: 1012019
    ## epsg (SRID):    32617
    ## proj4string:    +proj=utm +zone=17 +datum=WGS84 +units=m +no_defs
    ## First 10 features:
    ##    id Precipitation Age.cat Geology  Habitat Stream EnvHet UTM.EW  UTM.NS
    ## 1   1          2530      c3      Tb OldSlope    Yes 0.6272 625754 1011569
    ## 2   2          2530      c3      Tb   OldLow    Yes 0.3936 625754 1011669
    ## 3   3          2530      c3      Tb   OldLow     No 0.0000 625754 1011769
    ## 4   4          2530      c3      Tb   OldLow     No 0.0000 625754 1011869
    ## 5   5          2530      c3      Tb OldSlope     No 0.4608 625754 1011969
    ## 6   6          2530      c3      Tb   OldLow     No 0.0768 625854 1011569
    ## 7   7          2530      c3      Tb   OldLow    Yes 0.3808 625854 1011669
    ## 8   8          2530      c3      Tb   OldLow    Yes 0.2112 625854 1011769
    ## 9   9          2530      c3      Tb   OldLow     No 0.0000 625854 1011869
    ## 10 10          2530      c3      Tb   OldLow     No 0.0000 625854 1011969
    ##    elevation     convex    slope    aspectEW   aspectNS
    ## 1   130.2525  -7.872500 6.694828 -0.89108252 -0.4538413
    ## 2   136.8100 -10.700000 5.086842 -0.21903766 -0.9757164
    ## 3   143.6775 -14.667500 3.104794  0.03051372 -0.9995343
    ## 4   147.0075 -16.757500 1.872813 -0.86414183 -0.5032483
    ## 5   144.3850 -12.485000 5.118725 -0.67148116  0.7410216
    ## 6   136.8750  -9.685000 2.945532 -0.86532324 -0.5012142
    ## 7   139.2450   0.250000 3.128125 -0.76696195 -0.6416926
    ## 8   142.6425   0.188125 2.831006  0.40623791 -0.9137673
    ## 9   145.8750   1.129687 2.875696  0.89684900 -0.4423368
    ## 10  146.8500 -11.940000 3.607031  0.97247911  0.2329901
    ##                          geometry
    ## 1  POLYGON ((625704 1011519, 6...
    ## 2  POLYGON ((625704 1011619, 6...
    ## 3  POLYGON ((625704 1011719, 6...
    ## 4  POLYGON ((625704 1011819, 6...
    ## 5  POLYGON ((625704 1011919, 6...
    ## 6  POLYGON ((625804 1011519, 6...
    ## 7  POLYGON ((625804 1011619, 6...
    ## 8  POLYGON ((625804 1011719, 6...
    ## 9  POLYGON ((625804 1011819, 6...
    ## 10 POLYGON ((625804 1011919, 6...

### También podemos usar

Requiere que se haya cargado ya la colección tidyverse

``` r
censo_apcyn_melic_saptc %>% tibble
```

    ## # A tibble: 18,426 x 30
    ##    treeID stemID tag   StemTag sp    quadrat    gx    gy MeasureID CensusID
    ##     <int>  <int> <chr> <chr>   <chr> <chr>   <dbl> <dbl>     <int>    <int>
    ##  1     35      1 0000… ""      tri2… 4921     992.  432.        12      171
    ##  2     41      4 0000…  <NA>   guar… 4921     991.  420.        18      171
    ##  3     47      1 0000… ""      tri2… 4919     988.  392.        27      171
    ##  4     75      1 0000… ""      tri2… 4915     994.  308.        47      171
    ##  5     90      1 0000… ""      pout… 4912     987   250.        58      171
    ##  6     99      1 0000… ""      tab2… 4912     982.  241.        62      171
    ##  7    108      1 0000… ""      tri2… 4911     997   220.        68      171
    ##  8    113      1 0000… ""      chr2… 4909     998.  195.        71      171
    ##  9    125      1 0001… ""      guar… 4908     990.  168.        78      171
    ## 10    135      1 0001… ""      chr2… 4906     996.  134.        84      171
    ## # … with 18,416 more rows, and 20 more variables: dbh <dbl>, pom <chr>,
    ## #   hom <dbl>, ExactDate <chr>, DFstatus <chr>, codes <chr>,
    ## #   nostems <dbl>, date <dbl>, status <chr>, agb <dbl>, Latin <chr>,
    ## #   Genus <chr>, Species <chr>, Family <chr>, speciesID <int>,
    ## #   authority <chr>, IDlevel <chr>, syn <chr>, subsp <chr>, quad1ha <dbl>

``` r
mc_apcyn_melic_saptc %>% tibble
```

    ## # A tibble: 50 x 16
    ##    `Aspidosperma s… `Cedrela odorat… `Chrysophyllum … `Chrysophyllum …
    ##               <dbl>            <dbl>            <dbl>            <dbl>
    ##  1                3                0               21                2
    ##  2                2                0               11                1
    ##  3                2                0               19                1
    ##  4                3                0               38                2
    ##  5                2                1               21                0
    ##  6                3                0               18                0
    ##  7                0                0                9                3
    ##  8                4                0               14                0
    ##  9                5                0               25                0
    ## 10                1                0               17                1
    ## # … with 40 more rows, and 12 more variables: `Guarea bullata` <dbl>,
    ## #   `Guarea grandifolia` <dbl>, `Guarea guidonia` <dbl>, `Lacmellea
    ## #   panamensis` <dbl>, `Pouteria fossicola` <dbl>, `Pouteria
    ## #   reticulata` <dbl>, `Pouteria stipitata` <dbl>, `Rauvolfia
    ## #   littoralis` <dbl>, `Tabernaemontana arborea` <dbl>, `Thevetia
    ## #   ahouai` <dbl>, `Trichilia pallida` <dbl>, `Trichilia
    ## #   tuberculata` <dbl>

### Lista de especies

``` r
sort(colnames(mc_apcyn_melic_saptc))
```

    ##  [1] "Aspidosperma spruceanum" "Cedrela odorata"        
    ##  [3] "Chrysophyllum argenteum" "Chrysophyllum cainito"  
    ##  [5] "Guarea bullata"          "Guarea grandifolia"     
    ##  [7] "Guarea guidonia"         "Lacmellea panamensis"   
    ##  [9] "Pouteria fossicola"      "Pouteria reticulata"    
    ## [11] "Pouteria stipitata"      "Rauvolfia littoralis"   
    ## [13] "Tabernaemontana arborea" "Thevetia ahouai"        
    ## [15] "Trichilia pallida"       "Trichilia tuberculata"

### Número de sitios, tanto en matriz de comunidad como en ambiental

Verifica que
coinciden

``` r
nrow(mc_apcyn_melic_saptc) #En la matriz de comunidad
```

    ## [1] 50

``` r
nrow(bci_env_grid) #En la matriz ambiental
```

    ## [1] 50

### Riqueza numérica de especies (usando matriz de comunidad) por quadrat

Nota: cargar paquete vegan arriba, en el área de
    paquetes

``` r
specnumber(mc_apcyn_melic_saptc)
```

    ##  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 
    ## 12 11 12 12 14 11 11 11 11 11 12 12 10 13 12 11 11 12 11 12 12 11 13 11 12 
    ## 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 
    ##  8 11 13 12 12 10 12 12 11 11 11 11 11 11 11 12 13 10 13 10 11 11 12 10 12

``` r
sort(specnumber(mc_apcyn_melic_saptc)) # Ordenados ascendentemente
```

    ## 26 13 31 43 45 49  2  6  7  8  9 10 16 17 19 22 24 27 34 35 36 37 38 39 40 
    ##  8 10 10 10 10 10 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 
    ## 46 47  1  3  4 11 12 15 18 20 21 25 29 30 32 33 41 48 50 14 23 28 42 44  5 
    ## 11 11 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 13 13 13 13 13 14

``` r
summary(specnumber(mc_apcyn_melic_saptc)) # Resumen estadístico
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    8.00   11.00   11.00   11.44   12.00   14.00

### Abundancia de especies por quadrat

``` r
sort(rowSums(mc_apcyn_melic_saptc))
```

    ##  46   9   5   6  19  10   3   8  17  24  12  22   4   7  29  47  11  14 
    ## 167 180 188 193 201 214 215 215 217 219 231 234 237 253 254 257 266 278 
    ##  21  42  13   1  28   2  15  37  16  36  38  41  26  18  27  32  20  43 
    ## 282 298 306 309 311 313 315 340 349 350 359 361 371 385 414 426 429 430 
    ##  31  30  23  50  35  25  39  34  40  45  33  44  48  49 
    ## 432 459 491 494 531 551 555 581 593 609 610 684 724 745

``` r
summary(rowSums(mc_apcyn_melic_saptc)) # Resumen estadístico
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   167.0   241.0   327.5   368.5   452.2   745.0

### Abundancia por especie

``` r
sort(colSums(mc_apcyn_melic_saptc))
```

    ##    Rauvolfia littoralis      Pouteria fossicola         Cedrela odorata 
    ##                       1                       3                      12 
    ##      Pouteria stipitata      Guarea grandifolia         Thevetia ahouai 
    ##                      60                      65                      84 
    ##    Lacmellea panamensis   Chrysophyllum cainito       Trichilia pallida 
    ##                     102                     171                     472 
    ## Aspidosperma spruceanum Chrysophyllum argenteum          Guarea bullata 
    ##                     473                     711                     725 
    ##     Pouteria reticulata Tabernaemontana arborea         Guarea guidonia 
    ##                    1084                    1732                    1889 
    ##   Trichilia tuberculata 
    ##                   10842

``` r
summary(colSums(mc_apcyn_melic_saptc)) # Resumen estadístico
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ##     1.00    63.75   321.50  1151.62   814.75 10842.00

### Riqueza numérica de toda la “comunidad”

``` r
specnumber(colSums(mc_apcyn_melic_saptc))
```

    ## [1] 16

### Abundancia de toda la comunidad

``` r
sum(colSums(mc_apcyn_melic_saptc))
```

    ## [1] 18426

### Una tabla para el manuscrito, es necesario asignarle nombre

Para esto, usaré la colección “tidyverse”

``` r
abun_sp <- censo_apcyn_melic_saptc %>%
  group_by(Latin) %>% 
  count() %>% 
  arrange(desc(n))
abun_sp
```

    ## # A tibble: 16 x 2
    ## # Groups:   Latin [16]
    ##    Latin                       n
    ##    <chr>                   <int>
    ##  1 Trichilia tuberculata   10842
    ##  2 Guarea guidonia          1889
    ##  3 Tabernaemontana arborea  1732
    ##  4 Pouteria reticulata      1084
    ##  5 Guarea bullata            725
    ##  6 Chrysophyllum argenteum   711
    ##  7 Aspidosperma spruceanum   473
    ##  8 Trichilia pallida         472
    ##  9 Chrysophyllum cainito     171
    ## 10 Lacmellea panamensis      102
    ## 11 Thevetia ahouai            84
    ## 12 Guarea grandifolia         65
    ## 13 Pouteria stipitata         60
    ## 14 Cedrela odorata            12
    ## 15 Pouteria fossicola          3
    ## 16 Rauvolfia littoralis        1

### Un gráfico para el manuscrito

Gráfico de mosaicos de la abundancia por especie por
cuadros

``` r
abun_sp_q <- crear_grafico_mosaico_de_mc(mc_apcyn_melic_saptc, tam_rotulo = 6)
abun_sp_q
```

![](aed_1_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->
