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
  max.order = 2,
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
    ##  Selected number of species: 0 
    ##  Number of species associated to 1 group: 0 
    ## 
    ##  List of species associated to each combination: 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(p_upgma_adj <- p.adjust(iva_upgma_k2$sign$p.value))
```

    ##  [1]    NA    NA    NA    NA    NA    NA    NA    NA 0.774    NA    NA
    ## [12] 1.000    NA    NA    NA    NA

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
    ## Aspidosperma spruceanum 0.7798149 0.3608012
    ## Cedrela odorata         0.1143075 0.0000000
    ## Chrysophyllum argenteum 0.5534121 0.7547011
    ## Chrysophyllum cainito   0.7234123 0.2050809
    ## Guarea bullata          0.6163325 0.7177455
    ## Guarea grandifolia      0.4835455 0.3862027
    ## Guarea guidonia         0.7428970 0.5359806
    ## Lacmellea panamensis    0.4620704 0.7181044
    ## Pouteria fossicola      0.0000000 0.0000000
    ## Pouteria reticulata     0.6421125 0.6777045
    ## Pouteria stipitata      0.3441090 0.4812265
    ## Rauvolfia littoralis    0.0000000 0.0000000
    ## Tabernaemontana arborea 0.5724008 0.7492792
    ## Thevetia ahouai         0.1905811 0.0000000
    ## Trichilia pallida       0.5875200 0.5471645
    ## Trichilia tuberculata   0.8475516 0.4468437
    ## 
    ## $upperCI
    ##                                 1         2
    ## Aspidosperma spruceanum 0.9185934 0.5907738
    ## Cedrela odorata         0.4612656 0.7006490
    ## Chrysophyllum argenteum 0.6559482 0.8326664
    ## Chrysophyllum cainito   0.9077181 0.6280305
    ## Guarea bullata          0.6962305 0.7871800
    ## Guarea grandifolia      0.7280252 0.7664855
    ## Guarea guidonia         0.8441442 0.6692231
    ## Lacmellea panamensis    0.6407155 0.8358169
    ## Pouteria fossicola      0.3123475 0.6490940
    ## Pouteria reticulata     0.7352528 0.7665336
    ## Pouteria stipitata      0.6094883 0.8609161
    ## Rauvolfia littoralis    0.2886751 0.0000000
    ## Tabernaemontana arborea 0.6620528 0.8197764
    ## Thevetia ahouai         0.5577734 0.6479563
    ## Trichilia pallida       0.8205361 0.7939992
    ## Trichilia tuberculata   0.8944652 0.5305145

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
    ##  Selected number of species: 4 
    ##  Number of species associated to 1 group: 0 
    ##  Number of species associated to 2 groups: 4 
    ## 
    ##  List of species associated to each combination: 
    ## 
    ##  Group 1+2  #sps.  2 
    ##                              A      B  stat p.value    
    ## Tabernaemontana arborea 0.8014 1.0000 0.895   0.001 ***
    ## Chrysophyllum argenteum 0.7526 1.0000 0.868   0.048 *  
    ## 
    ##  Group 1+3  #sps.  1 
    ##                            A      B  stat p.value    
    ## Trichilia tuberculata 0.8602 1.0000 0.927   0.001 ***
    ## 
    ##  Group 2+3  #sps.  1 
    ##                     A      B  stat p.value  
    ## Guarea bullata 0.7552 1.0000 0.869   0.039 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(p_ward_adj <- p.adjust(iva_ward_k3$sign$p.value))
```

    ##  [1] 1.000 1.000 0.624 1.000 0.546 1.000 1.000 1.000 1.000 1.000 1.000
    ## [12] 1.000 0.016 1.000 1.000 0.016

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
    ## Aspidosperma spruceanum 0.3370513 0.2904383 0.70531410
    ## Cedrela odorata         0.0000000 0.0000000 0.09022199
    ## Chrysophyllum argenteum 0.4712062 0.6345897 0.44539996
    ## Chrysophyllum cainito   0.4279181 0.1632022 0.53435079
    ## Guarea bullata          0.4456249 0.6003775 0.54339433
    ## Guarea grandifolia      0.3936750 0.4986169 0.25254790
    ## Guarea guidonia         0.4779667 0.4420826 0.61314514
    ## Lacmellea panamensis    0.3796459 0.6259240 0.33049591
    ## Pouteria fossicola      0.0000000 0.0000000 0.00000000
    ## Pouteria reticulata     0.5197860 0.5311765 0.51110282
    ## Pouteria stipitata      0.3459600 0.2363597 0.21182964
    ## Rauvolfia littoralis    0.0000000 0.0000000 0.00000000
    ## Tabernaemontana arborea 0.5333864 0.6234150 0.41140254
    ## Thevetia ahouai         0.1376653 0.0000000 0.02194744
    ## Trichilia pallida       0.4158020 0.4837578 0.46385645
    ## Trichilia tuberculata   0.4659942 0.3323085 0.74561335
    ## 
    ## $upperCI
    ##                                 1         2         3
    ## Aspidosperma spruceanum 0.5701113 0.5007403 0.8585663
    ## Cedrela odorata         0.2806849 0.8192597 0.5002127
    ## Chrysophyllum argenteum 0.5803989 0.7438092 0.5491941
    ## Chrysophyllum cainito   0.6721638 0.5553108 0.7351052
    ## Guarea bullata          0.5478293 0.6856766 0.6230853
    ## Guarea grandifolia      0.6788829 0.6798425 0.5207770
    ## Guarea guidonia         0.5982204 0.5750475 0.7340778
    ## Lacmellea panamensis    0.6041190 0.7461361 0.5181044
    ## Pouteria fossicola      0.4588315 0.7453560 0.0000000
    ## Pouteria reticulata     0.6217447 0.6640449 0.6053032
    ## Pouteria stipitata      0.6620830 0.7709951 0.5031314
    ## Rauvolfia littoralis    0.0000000 0.0000000 0.3692745
    ## Tabernaemontana arborea 0.6382067 0.7201566 0.4795153
    ## Thevetia ahouai         0.6179144 0.7040077 0.3478466
    ## Trichilia pallida       0.6407842 0.7226901 0.6553429
    ## Trichilia tuberculata   0.5464102 0.4093703 0.8059659

## Análisis de especies con preferencia por hábitat mediante el coeficiente de correlación biserial puntual

### UPGMA

``` r
phi_upgma_k2 <- multipatt(
  mi_fam,
  grupos_upgma_k2,
  func = "r.g",
  max.order = 2,
  control = how(nperm = 999))
summary(phi_upgma_k2, alpha = 0.01)
```

    ## 
    ##  Multilevel pattern analysis
    ##  ---------------------------
    ## 
    ##  Association function: r.g
    ##  Significance level (alpha): 0.01
    ## 
    ##  Total number of species: 16
    ##  Selected number of species: 3 
    ##  Number of species associated to 1 group: 3 
    ## 
    ##  List of species associated to each combination: 
    ## 
    ##  Group 1  #sps.  1 
    ##                        stat p.value    
    ## Trichilia tuberculata 0.676   0.001 ***
    ## 
    ##  Group 2  #sps.  2 
    ##                          stat p.value    
    ## Chrysophyllum argenteum 0.579   0.001 ***
    ## Tabernaemontana arborea 0.562   0.002 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

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
    ## Aspidosperma spruceanum  0.2956753 -0.57545640
    ## Cedrela odorata         -0.4731700 -0.29397237
    ## Chrysophyllum argenteum -0.7388810  0.43996221
    ## Chrysophyllum cainito    0.2459825 -0.64350120
    ## Guarea bullata          -0.6690382  0.07563813
    ## Guarea grandifolia      -0.2675323 -0.34783780
    ## Guarea guidonia          0.1963678 -0.59296114
    ## Lacmellea panamensis    -0.6647085  0.06270849
    ## Pouteria fossicola      -0.5000000 -0.22941573
    ## Pouteria reticulata     -0.5381022 -0.30838288
    ## Pouteria stipitata      -0.5803131 -0.19425717
    ## Rauvolfia littoralis     0.0000000 -0.21320072
    ## Tabernaemontana arborea -0.7862856  0.32240636
    ## Thevetia ahouai         -0.3005007 -0.28734759
    ## Trichilia pallida       -0.3809710 -0.46279985
    ## Trichilia tuberculata    0.6037522 -0.75559021
    ## 
    ## $upperCI
    ##                                   1          2
    ## Aspidosperma spruceanum  0.57504572 -0.2978240
    ## Cedrela odorata          0.29397237  0.4713132
    ## Chrysophyllum argenteum -0.44167394  0.7387060
    ## Chrysophyllum cainito    0.64215338 -0.2492651
    ## Guarea bullata          -0.07936946  0.6684915
    ## Guarea grandifolia       0.34689260  0.2641644
    ## Guarea guidonia          0.59268602 -0.1989227
    ## Lacmellea panamensis    -0.06430684  0.6646446
    ## Pouteria fossicola       0.22645541  0.5000000
    ## Pouteria reticulata      0.30037545  0.5367369
    ## Pouteria stipitata       0.19329920  0.5728842
    ## Rauvolfia littoralis     0.20555661  0.0000000
    ## Tabernaemontana arborea -0.32400469  0.7847987
    ## Thevetia ahouai          0.28712177  0.2933258
    ## Trichilia pallida        0.46254676  0.3736349
    ## Trichilia tuberculata    0.75468287 -0.6047774

Ward

``` r
phi_ward_k3 <- multipatt(
  mi_fam,
  grupos_ward_k3,
  func = "r.g",
  max.order = 2,
  control = how(nperm = 999))
summary(phi_ward_k3, alpha = 0.05)
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
    ## Aspidosperma spruceanum 0.520   0.005 ** 
    ## Guarea guidonia         0.414   0.040 *  
    ## 
    ##  Group 1+2  #sps.  1 
    ##                         stat p.value   
    ## Tabernaemontana arborea 0.61   0.002 **
    ## 
    ##  Group 2+3  #sps.  1 
    ##                 stat p.value  
    ## Guarea bullata 0.475   0.036 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

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
    ## Aspidosperma spruceanum -0.3810717 -0.4658156  0.342474353
    ## Cedrela odorata         -0.4370934 -0.2413899 -0.324604299
    ## Chrysophyllum argenteum -0.4333855  0.3953495 -0.534073498
    ## Chrysophyllum cainito   -0.1989332 -0.5496328 -0.006548647
    ## Guarea bullata          -0.6663942  0.1713286 -0.216759929
    ## Guarea grandifolia      -0.1574209 -0.2110079 -0.442291490
    ## Guarea guidonia         -0.3514338 -0.4529795  0.181502085
    ## Lacmellea panamensis    -0.3649330  0.2222787 -0.567108820
    ## Pouteria fossicola      -0.3123790 -0.1889822 -0.418330013
    ## Pouteria reticulata     -0.3549829 -0.3572597 -0.411309276
    ## Pouteria stipitata      -0.2526181 -0.3688541 -0.448413977
    ## Rauvolfia littoralis    -0.1524986 -0.1524986  0.000000000
    ## Tabernaemontana arborea -0.2011049  0.2397294 -0.713333312
    ## Thevetia ahouai         -0.2024183 -0.2325441 -0.313677107
    ## Trichilia pallida       -0.4047980 -0.2993815 -0.291684718
    ## Trichilia tuberculata   -0.3205511 -0.6453552  0.714880958
    ## 
    ## $upperCI
    ##                                    1           2           3
    ## Aspidosperma spruceanum  0.009655634 -0.15566325  0.67187485
    ## Cedrela odorata          0.061765775  0.66712438  0.33896557
    ## Chrysophyllum argenteum  0.047254082  0.79350808 -0.18715570
    ## Chrysophyllum cainito    0.307353562 -0.09193373  0.50170733
    ## Guarea bullata          -0.235344445  0.66461099  0.29403418
    ## Guarea grandifolia       0.418811522  0.33086712  0.08313122
    ## Guarea guidonia          0.081515313 -0.01123553  0.59221315
    ## Lacmellea panamensis     0.221878440  0.66213537 -0.12680586
    ## Pouteria fossicola       0.377964473  0.68041382  0.00000000
    ## Pouteria reticulata      0.366848097  0.64839864  0.24411753
    ## Pouteria stipitata       0.399470176  0.55939055  0.15616068
    ## Rauvolfia littoralis     0.000000000  0.00000000  0.30151134
    ## Tabernaemontana arborea  0.384416389  0.78839277 -0.49375568
    ## Thevetia ahouai          0.364259480  0.38644464  0.15965291
    ## Trichilia pallida        0.252103500  0.53589683  0.32912767
    ## Trichilia tuberculata   -0.112544888 -0.50592280  0.86614349
