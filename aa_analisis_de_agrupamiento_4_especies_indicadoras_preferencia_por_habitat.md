Análisis de agrupamiento (cluster analysis). <br> Parte 4: Especies
indicadoras, especies con preferencia por hábitats
================
JR
15 de noviembre, 2020

``` r
knitr::opts_chunk$set(fig.width=12, fig.height=8)
```

## Preámbulo

### Cargar paquetes

``` r
library(indicspecies)
```

    ## Loading required package: permute

``` r
source('biodata/funciones.R')
```

### Cargar datos

``` r
load('biodata/Apocynaceae-Meliaceae-Sapotaceae.Rdata')
mi_fam <- mc_apcyn_melic_saptc
grupos_upgma_k2 <- readRDS('grupos_upgma_k2.RDS')
table(grupos_upgma_k2)
```

    ## grupos_upgma_k2
    ##  1  2 
    ## 43  7

``` r
grupos_ward_k3 <- readRDS('grupos_ward_k3.RDS')
table(grupos_ward_k3)
```

    ## grupos_ward_k3
    ##  1  2  3 
    ## 20  5 25

## Análisis de especies indicadoras mediante IndVal

### UPGMA

``` r
iva_upgma_k2 <- multipatt(
  x = mi_fam,
  cluster = grupos_upgma_k2,
  func = 'IndVal.g',
  max.order = 1,
  control = how(nperm = 999))
summary(iva_upgma_k2, indvalcomp = TRUE)
```

    ## 
    ##  Multilevel pattern analysis
    ##  ---------------------------
    ## 
    ##  Association function: IndVal.g
    ##  Significance level (alpha): 0.05
    ## 
    ##  Total number of species: 16
    ##  Selected number of species: 7 
    ##  Number of species associated to 1 group: 7 
    ## 
    ##  List of species associated to each combination: 
    ## 
    ##  Group 1  #sps.  4 
    ##                              A      B  stat p.value    
    ## Trichilia tuberculata   0.7629 1.0000 0.873   0.001 ***
    ## Aspidosperma spruceanum 0.7780 0.9535 0.861   0.016 *  
    ## Chrysophyllum cainito   0.7238 0.9070 0.810   0.032 *  
    ## Guarea guidonia         0.6344 1.0000 0.797   0.033 *  
    ## 
    ##  Group 2  #sps.  3 
    ##                              A      B  stat p.value  
    ## Chrysophyllum argenteum 0.6294 1.0000 0.793   0.011 *
    ## Tabernaemontana arborea 0.6196 1.0000 0.787   0.021 *
    ## Lacmellea panamensis    0.6143 1.0000 0.784   0.044 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
colSums(mi_fam)
```

    ## Aspidosperma spruceanum         Cedrela odorata Chrysophyllum argenteum 
    ##                     473                      12                     711 
    ##   Chrysophyllum cainito          Guarea bullata      Guarea grandifolia 
    ##                     171                     725                      65 
    ##         Guarea guidonia    Lacmellea panamensis      Pouteria fossicola 
    ##                    1889                     102                       3 
    ##     Pouteria reticulata      Pouteria stipitata    Rauvolfia littoralis 
    ##                    1084                      60                       1 
    ## Tabernaemontana arborea         Thevetia ahouai       Trichilia pallida 
    ##                    1732                      84                     472 
    ##   Trichilia tuberculata 
    ##                   10842

``` r
(p_upgma_adj <- p.adjust(iva_upgma_k2$sign$p.value))
```

    ##  [1] 0.224 1.000 0.165 0.384 1.000 1.000 0.384 0.440 1.000 1.000 1.000
    ## [12] 1.000 0.273 1.000 1.000 0.016

``` r
(iva_upgma_boot <- strassoc(
  X = mi_fam,
  cluster = grupos_upgma_k2,
  func = "IndVal.g",
  nboot = 1000))
```

    ## $stat
    ##                                 1         2
    ## Aspidosperma spruceanum 0.8612691 0.4712016
    ## Cedrela odorata         0.2702723 0.3968742
    ## Chrysophyllum argenteum 0.6087628 0.7933523
    ## Chrysophyllum cainito   0.8102438 0.4441462
    ## Guarea bullata          0.6551811 0.7554719
    ## Guarea grandifolia      0.6018064 0.6299577
    ## Guarea guidonia         0.7965090 0.6046266
    ## Lacmellea panamensis    0.5522530 0.7837638
    ## Pouteria fossicola      0.1068827 0.3282825
    ## Pouteria reticulata     0.6882979 0.7254281
    ## Pouteria stipitata      0.4788347 0.7204961
    ## Rauvolfia littoralis    0.1524986 0.0000000
    ## Tabernaemontana arborea 0.6167646 0.7871477
    ## Thevetia ahouai         0.3862963 0.3199702
    ## Trichilia pallida       0.7030733 0.6939562
    ## Trichilia tuberculata   0.8734660 0.4868852
    ## 
    ## $lowerCI
    ##                                 1         2
    ## Aspidosperma spruceanum 0.7813540 0.3671940
    ## Cedrela odorata         0.1147638 0.0000000
    ## Chrysophyllum argenteum 0.5495272 0.7536406
    ## Chrysophyllum cainito   0.7083364 0.2084486
    ## Guarea bullata          0.6127778 0.7193906
    ## Guarea grandifolia      0.4827641 0.4092753
    ## Guarea guidonia         0.7410722 0.5435370
    ## Lacmellea panamensis    0.4595519 0.7185432
    ## Pouteria fossicola      0.0000000 0.0000000
    ## Pouteria reticulata     0.6392701 0.6730598
    ## Pouteria stipitata      0.3416718 0.4847416
    ## Rauvolfia littoralis    0.0000000 0.0000000
    ## Tabernaemontana arborea 0.5702817 0.7481751
    ## Thevetia ahouai         0.2042879 0.0000000
    ## Trichilia pallida       0.5815853 0.5456596
    ## Trichilia tuberculata   0.8480616 0.4457490
    ## 
    ## $upperCI
    ##                                 1         2
    ## Aspidosperma spruceanum 0.9184514 0.5889780
    ## Cedrela odorata         0.4662524 0.7146415
    ## Chrysophyllum argenteum 0.6572711 0.8347446
    ## Chrysophyllum cainito   0.8991574 0.6255432
    ## Guarea bullata          0.6944361 0.7899383
    ## Guarea grandifolia      0.7124601 0.7638765
    ## Guarea guidonia         0.8393480 0.6708204
    ## Lacmellea panamensis    0.6447197 0.8355894
    ## Pouteria fossicola      0.3261640 0.6650622
    ## Pouteria reticulata     0.7392821 0.7688097
    ## Pouteria stipitata      0.6178153 0.8578164
    ## Rauvolfia littoralis    0.2738613 0.0000000
    ## Tabernaemontana arborea 0.6633678 0.8212148
    ## Thevetia ahouai         0.5657995 0.6435382
    ## Trichilia pallida       0.8148217 0.7982979
    ## Trichilia tuberculata   0.8950264 0.5298364

Ward

``` r
iva_ward_k3 <- multipatt(
  x = mi_fam,
  cluster = grupos_ward_k3,
  func = 'IndVal.g',
  max.order = 2,
  control = how(nperm = 999))
summary(iva_ward_k3, indvalcomp = TRUE)
```

    ## 
    ##  Multilevel pattern analysis
    ##  ---------------------------
    ## 
    ##  Association function: IndVal.g
    ##  Significance level (alpha): 0.05
    ## 
    ##  Total number of species: 16
    ##  Selected number of species: 2 
    ##  Number of species associated to 1 group: 0 
    ##  Number of species associated to 2 groups: 2 
    ## 
    ##  List of species associated to each combination: 
    ## 
    ##  Group 1+2  #sps.  1 
    ##                              A      B  stat p.value   
    ## Tabernaemontana arborea 0.8014 1.0000 0.895   0.006 **
    ## 
    ##  Group 1+3  #sps.  1 
    ##                            A      B  stat p.value   
    ## Trichilia tuberculata 0.8602 1.0000 0.927   0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
colSums(mi_fam)
```

    ## Aspidosperma spruceanum         Cedrela odorata Chrysophyllum argenteum 
    ##                     473                      12                     711 
    ##   Chrysophyllum cainito          Guarea bullata      Guarea grandifolia 
    ##                     171                     725                      65 
    ##         Guarea guidonia    Lacmellea panamensis      Pouteria fossicola 
    ##                    1889                     102                       3 
    ##     Pouteria reticulata      Pouteria stipitata    Rauvolfia littoralis 
    ##                    1084                      60                       1 
    ## Tabernaemontana arborea         Thevetia ahouai       Trichilia pallida 
    ##                    1732                      84                     472 
    ##   Trichilia tuberculata 
    ##                   10842

``` r
(p_ward_adj <- p.adjust(iva_ward_k3$sign$p.value))
```

    ##  [1] 0.948 1.000 0.756 1.000 0.756 1.000 1.000 1.000 1.000 1.000 1.000
    ## [12] 1.000 0.090 1.000 1.000 0.032

``` r
(iva_ward_boot <- strassoc(
  X = mi_fam,
  cluster = grupos_ward_k3,
  func = "IndVal.g",
  nboot = 1000))
```

    ## $stat
    ##                                 1         2         3
    ## Aspidosperma spruceanum 0.4542568 0.3871845 0.7879158
    ## Cedrela odorata         0.1104315 0.4417261 0.2793721
    ## Chrysophyllum argenteum 0.5285637 0.6879354 0.4973583
    ## Chrysophyllum cainito   0.5522071 0.4017897 0.6358498
    ## Guarea bullata          0.4947575 0.6459093 0.5813916
    ## Guarea grandifolia      0.5450043 0.5886719 0.3930082
    ## Guarea guidonia         0.5375324 0.5042113 0.6758920
    ## Lacmellea panamensis    0.4930456 0.6852833 0.4256941
    ## Pouteria fossicola      0.1825742 0.3651484 0.0000000
    ## Pouteria reticulata     0.5693960 0.6023946 0.5593827
    ## Pouteria stipitata      0.4974683 0.5685352 0.3550501
    ## Rauvolfia littoralis    0.0000000 0.0000000 0.2000000
    ## Tabernaemontana arborea 0.5886238 0.6744480 0.4456926
    ## Thevetia ahouai         0.4082483 0.3333333 0.1632993
    ## Trichilia pallida       0.5314437 0.6151247 0.5579887
    ## Trichilia tuberculata   0.5049846 0.3738488 0.7779638
    ## 
    ## $lowerCI
    ##                                 1         2          3
    ## Aspidosperma spruceanum 0.3346646 0.2863689 0.70759984
    ## Cedrela odorata         0.0000000 0.0000000 0.08987182
    ## Chrysophyllum argenteum 0.4782818 0.6361843 0.44706410
    ## Chrysophyllum cainito   0.4190203 0.1675416 0.53078379
    ## Guarea bullata          0.4474066 0.5996285 0.54090219
    ## Guarea grandifolia      0.3996313 0.4992872 0.26087943
    ## Guarea guidonia         0.4776207 0.4422923 0.62051911
    ## Lacmellea panamensis    0.3744990 0.6350803 0.32915005
    ## Pouteria fossicola      0.0000000 0.0000000 0.00000000
    ## Pouteria reticulata     0.5206010 0.5314583 0.51468277
    ## Pouteria stipitata      0.3405835 0.2711718 0.21123850
    ## Rauvolfia littoralis    0.0000000 0.0000000 0.00000000
    ## Tabernaemontana arborea 0.5369224 0.6247982 0.41314051
    ## Thevetia ahouai         0.1226791 0.0000000 0.01685853
    ## Trichilia pallida       0.4215325 0.4812026 0.45822963
    ## Trichilia tuberculata   0.4691351 0.3294451 0.74572488
    ## 
    ## $upperCI
    ##                                 1         2         3
    ## Aspidosperma spruceanum 0.5611074 0.4948717 0.8572863
    ## Cedrela odorata         0.2881952 0.8289085 0.5051435
    ## Chrysophyllum argenteum 0.5803257 0.7409525 0.5431567
    ## Chrysophyllum cainito   0.6645650 0.5565966 0.7274623
    ## Guarea bullata          0.5424890 0.6839698 0.6216221
    ## Guarea grandifolia      0.6764814 0.6762352 0.5260766
    ## Guarea guidonia         0.5940423 0.5664224 0.7323037
    ## Lacmellea panamensis    0.5928105 0.7457545 0.5228315
    ## Pouteria fossicola      0.4662524 0.7377111 0.0000000
    ## Pouteria reticulata     0.6189614 0.6656997 0.6059755
    ## Pouteria stipitata      0.6508483 0.7673233 0.5096472
    ## Rauvolfia littoralis    0.0000000 0.0000000 0.3692745
    ## Tabernaemontana arborea 0.6361285 0.7184562 0.4824473
    ## Thevetia ahouai         0.6154060 0.6892173 0.3567758
    ## Trichilia pallida       0.6454972 0.7226310 0.6563249
    ## Trichilia tuberculata   0.5429037 0.4095112 0.8075115

## Análisis de especies con preferencia por hábitat mediante el coeficiente de correlación biserial puntual

### UPGMA

``` r
phi_upgma_k2 <- multipatt(
  mi_fam,
  grupos_upgma_k2,
  func = "r.g",
  max.order = 1,
  control = how(nperm = 999))
summary(phi_upgma_k2)
```

    ## 
    ##  Multilevel pattern analysis
    ##  ---------------------------
    ## 
    ##  Association function: r.g
    ##  Significance level (alpha): 0.05
    ## 
    ##  Total number of species: 16
    ##  Selected number of species: 6 
    ##  Number of species associated to 1 group: 6 
    ## 
    ##  List of species associated to each combination: 
    ## 
    ##  Group 1  #sps.  4 
    ##                          stat p.value    
    ## Trichilia tuberculata   0.676   0.001 ***
    ## Chrysophyllum cainito   0.457   0.022 *  
    ## Aspidosperma spruceanum 0.446   0.022 *  
    ## Guarea guidonia         0.409   0.032 *  
    ## 
    ##  Group 2  #sps.  2 
    ##                          stat p.value    
    ## Chrysophyllum argenteum 0.579   0.002 ** 
    ## Tabernaemontana arborea 0.562   0.001 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
colSums(mi_fam)
```

    ## Aspidosperma spruceanum         Cedrela odorata Chrysophyllum argenteum 
    ##                     473                      12                     711 
    ##   Chrysophyllum cainito          Guarea bullata      Guarea grandifolia 
    ##                     171                     725                      65 
    ##         Guarea guidonia    Lacmellea panamensis      Pouteria fossicola 
    ##                    1889                     102                       3 
    ##     Pouteria reticulata      Pouteria stipitata    Rauvolfia littoralis 
    ##                    1084                      60                       1 
    ## Tabernaemontana arborea         Thevetia ahouai       Trichilia pallida 
    ##                    1732                      84                     472 
    ##   Trichilia tuberculata 
    ##                   10842

``` r
(phi_upgma_boot <- strassoc(
  X = mi_fam,
  cluster = grupos_upgma_k2,
  func = "r.g",
  nboot = 1000))
```

    ## $stat
    ##                                   1           2
    ## Aspidosperma spruceanum  0.44622874 -0.44622874
    ## Cedrela odorata         -0.04616033  0.04616033
    ## Chrysophyllum argenteum -0.57877980  0.57877980
    ## Chrysophyllum cainito    0.45652692 -0.45652692
    ## Guarea bullata          -0.39527724  0.39527724
    ## Guarea grandifolia       0.08088926 -0.08088926
    ## Guarea guidonia          0.40927828 -0.40927828
    ## Lacmellea panamensis    -0.38757239  0.38757239
    ## Pouteria fossicola      -0.16453652  0.16453652
    ## Pouteria reticulata     -0.15119392  0.15119392
    ## Pouteria stipitata      -0.22198686  0.22198686
    ## Rauvolfia littoralis     0.10846523 -0.10846523
    ## Tabernaemontana arborea -0.56156294  0.56156294
    ## Thevetia ahouai          0.09277495 -0.09277495
    ## Trichilia pallida        0.04645371 -0.04645371
    ## Trichilia tuberculata    0.67609388 -0.67609388
    ## 
    ## $lowerCI
    ##                                  1           2
    ## Aspidosperma spruceanum  0.2852706 -0.57180700
    ## Cedrela odorata         -0.4883149 -0.28697202
    ## Chrysophyllum argenteum -0.7487606  0.43864385
    ## Chrysophyllum cainito    0.2526679 -0.63988734
    ## Guarea bullata          -0.6502324  0.12203381
    ## Guarea grandifolia      -0.2284558 -0.34982151
    ## Guarea guidonia          0.1996681 -0.58397902
    ## Lacmellea panamensis    -0.6640044  0.04788912
    ## Pouteria fossicola      -0.4628529 -0.22941573
    ## Pouteria reticulata     -0.5315693 -0.26468319
    ## Pouteria stipitata      -0.5868994 -0.16920257
    ## Rauvolfia littoralis     0.0000000 -0.19738551
    ## Tabernaemontana arborea -0.7798319  0.32106440
    ## Thevetia ahouai         -0.3074210 -0.28508085
    ## Trichilia pallida       -0.3524051 -0.46896206
    ## Trichilia tuberculata    0.6017985 -0.74792709
    ## 
    ## $upperCI
    ##                                   1          2
    ## Aspidosperma spruceanum  0.57178575 -0.2865875
    ## Cedrela odorata          0.28697202  0.4883149
    ## Chrysophyllum argenteum -0.43992839  0.7481080
    ## Chrysophyllum cainito    0.63985587 -0.2566526
    ## Guarea bullata          -0.12271671  0.6497919
    ## Guarea grandifolia       0.34530507  0.2259852
    ## Guarea guidonia          0.58054051 -0.2017531
    ## Lacmellea panamensis    -0.04858339  0.6609383
    ## Pouteria fossicola       0.22941573  0.4628529
    ## Pouteria reticulata      0.25725353  0.5270906
    ## Pouteria stipitata       0.16807415  0.5841769
    ## Rauvolfia littoralis     0.19738551  0.0000000
    ## Tabernaemontana arborea -0.32122790  0.7779341
    ## Thevetia ahouai          0.28485729  0.3019397
    ## Trichilia pallida        0.46802565  0.3449822
    ## Trichilia tuberculata    0.74784894 -0.6057331

Ward

``` r
phi_ward_k3 <- multipatt(
  mi_fam,
  grupos_ward_k3,
  func = "r.g",
  max.order = 2,
  control = how(nperm = 999))
summary(phi_ward_k3)
```

    ## 
    ##  Multilevel pattern analysis
    ##  ---------------------------
    ## 
    ##  Association function: r.g
    ##  Significance level (alpha): 0.05
    ## 
    ##  Total number of species: 16
    ##  Selected number of species: 7 
    ##  Number of species associated to 1 group: 5 
    ##  Number of species associated to 2 groups: 2 
    ## 
    ##  List of species associated to each combination: 
    ## 
    ##  Group 2  #sps.  2 
    ##                          stat p.value   
    ## Chrysophyllum argenteum 0.613   0.002 **
    ## Lacmellea panamensis    0.460   0.038 * 
    ## 
    ##  Group 3  #sps.  3 
    ##                          stat p.value    
    ## Trichilia tuberculata   0.794   0.001 ***
    ## Aspidosperma spruceanum 0.520   0.007 ** 
    ## Guarea guidonia         0.414   0.031 *  
    ## 
    ##  Group 1+2  #sps.  1 
    ##                         stat p.value   
    ## Tabernaemontana arborea 0.61   0.005 **
    ## 
    ##  Group 2+3  #sps.  1 
    ##                 stat p.value  
    ## Guarea bullata 0.475   0.029 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
colSums(mi_fam)
```

    ## Aspidosperma spruceanum         Cedrela odorata Chrysophyllum argenteum 
    ##                     473                      12                     711 
    ##   Chrysophyllum cainito          Guarea bullata      Guarea grandifolia 
    ##                     171                     725                      65 
    ##         Guarea guidonia    Lacmellea panamensis      Pouteria fossicola 
    ##                    1889                     102                       3 
    ##     Pouteria reticulata      Pouteria stipitata    Rauvolfia littoralis 
    ##                    1084                      60                       1 
    ## Tabernaemontana arborea         Thevetia ahouai       Trichilia pallida 
    ##                    1732                      84                     472 
    ##   Trichilia tuberculata 
    ##                   10842

``` r
(phi_ward_boot <- strassoc(
  X = mi_fam,
  cluster = grupos_ward_k3,
  func = "r.g",
  nboot = 1000))
```

    ## $stat
    ##                                   1           2           3
    ## Aspidosperma spruceanum -0.18806952 -0.33151238  0.51958190
    ## Cedrela odorata         -0.20466732  0.14956458  0.05510274
    ## Chrysophyllum argenteum -0.23639113  0.61304873 -0.37665761
    ## Chrysophyllum cainito    0.06064008 -0.31390158  0.25326150
    ## Guarea bullata          -0.47474811  0.44964124  0.02510687
    ## Guarea grandifolia       0.14725585  0.03100123 -0.17825709
    ## Guarea guidonia         -0.14895761 -0.26543327  0.41439088
    ## Lacmellea panamensis    -0.09955764  0.46045409 -0.36089645
    ## Pouteria fossicola       0.00000000  0.23570226 -0.23570226
    ## Pouteria reticulata     -0.05272131  0.17077120 -0.11804989
    ## Pouteria stipitata       0.03975952  0.13915832 -0.17891784
    ## Rauvolfia littoralis    -0.08219949 -0.08219949  0.16439899
    ## Tabernaemontana arborea  0.05951128  0.55029413 -0.60980541
    ## Thevetia ahouai          0.16066502 -0.04016626 -0.12049877
    ## Trichilia pallida       -0.10141828  0.12677286 -0.02535457
    ## Trichilia tuberculata   -0.22883883 -0.56555453  0.79439336
    ## 
    ## $lowerCI
    ##                                  1          2            3
    ## Aspidosperma spruceanum -0.3690424 -0.4534220  0.332181218
    ## Cedrela odorata         -0.4406233 -0.2212524 -0.301735409
    ## Chrysophyllum argenteum -0.4376245  0.4345974 -0.559994784
    ## Chrysophyllum cainito   -0.2112656 -0.5206719 -0.008829348
    ## Guarea bullata          -0.6656606  0.1862935 -0.227964376
    ## Guarea grandifolia      -0.1307664 -0.2083868 -0.454603500
    ## Guarea guidonia         -0.3600590 -0.4351491  0.189246158
    ## Lacmellea panamensis    -0.3794816  0.2428361 -0.560027496
    ## Pouteria fossicola      -0.3162278 -0.1889822 -0.411113226
    ## Pouteria reticulata     -0.3480698 -0.3171528 -0.406773656
    ## Pouteria stipitata      -0.2687544 -0.3180328 -0.454907693
    ## Rauvolfia littoralis    -0.1507557 -0.1507557  0.000000000
    ## Tabernaemontana arborea -0.1876511  0.2984195 -0.715193747
    ## Thevetia ahouai         -0.2023950 -0.2298729 -0.335530206
    ## Trichilia pallida       -0.4152433 -0.2487789 -0.326354268
    ## Trichilia tuberculata   -0.3214949 -0.6374619  0.712807147
    ## 
    ## $upperCI
    ##                                   1           2           3
    ## Aspidosperma spruceanum  0.01444028 -0.15347342  0.67702625
    ## Cedrela odorata          0.01972596  0.63733435  0.33404844
    ## Chrysophyllum argenteum  0.01689096  0.80955625 -0.20629141
    ## Chrysophyllum cainito    0.32841007 -0.09357760  0.49709377
    ## Guarea bullata          -0.23978797  0.67083643  0.26785524
    ## Guarea grandifolia       0.42358904  0.30753179  0.08164415
    ## Guarea guidonia          0.08196809 -0.01991513  0.58170727
    ## Lacmellea panamensis     0.20696757  0.68843635 -0.12680872
    ## Pouteria fossicola       0.37796447  0.63245553  0.00000000
    ## Pouteria reticulata      0.29854327  0.63183654  0.20185150
    ## Pouteria stipitata       0.39098735  0.60192775  0.14282101
    ## Rauvolfia littoralis     0.00000000  0.00000000  0.30151134
    ## Tabernaemontana arborea  0.32696357  0.78922613 -0.51552611
    ## Thevetia ahouai          0.38200081  0.38589023  0.15326788
    ## Trichilia pallida        0.21590137  0.55900437  0.31792226
    ## Trichilia tuberculata   -0.11244564 -0.50468007  0.86768054
