Análisis espacial de datos ecológicos <br> Autocorrelación
================
JR
5 de diciembre, 2020

``` r
knitr::opts_chunk$set(fig.width=12, fig.height=8)
```

## Preámbulo

### Cargar paquetes

``` r
library(ape)
library(spdep)
```

    ## Loading required package: sp

    ## Loading required package: spData

    ## To access larger datasets in this package, install the spDataLarge
    ## package with: `install.packages('spDataLarge',
    ## repos='https://nowosad.github.io/drat/', type='source')`

    ## Loading required package: sf

    ## Linking to GEOS 3.6.2, GDAL 2.2.3, PROJ 4.9.3

    ## Registered S3 method overwritten by 'spdep':
    ##   method   from
    ##   plot.mst ape

``` r
library(ade4)
```

    ## 
    ## Attaching package: 'ade4'

    ## The following object is masked from 'package:spdep':
    ## 
    ##     mstree

``` r
library(adegraphics)
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

    ## 
    ## Attaching package: 'adegraphics'

    ## The following objects are masked from 'package:ade4':
    ## 
    ##     kplotsepan.coa, s.arrow, s.class, s.corcircle, s.distri,
    ##     s.image, s.label, s.logo, s.match, s.traject, s.value,
    ##     table.value, triangle.class

    ## The following object is masked from 'package:ape':
    ## 
    ##     zoom

``` r
library(adespatial)
```

    ## Registered S3 methods overwritten by 'adespatial':
    ##   method             from       
    ##   plot.multispati    adegraphics
    ##   print.multispati   ade4       
    ##   summary.multispati ade4

    ## 
    ## Attaching package: 'adespatial'

    ## The following object is masked from 'package:ade4':
    ## 
    ##     multispati

``` r
library(vegan)
```

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.5-6

``` r
library(tidyverse)
```

    ## ── Attaching packages ──────────────────────────────── tidyverse 1.2.1 ──

    ## ✓ ggplot2 3.3.2     ✓ purrr   0.3.4
    ## ✓ tibble  3.0.3     ✓ dplyr   0.8.3
    ## ✓ tidyr   1.0.0     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.4.0

    ## ── Conflicts ─────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(sf)
library(gridExtra)
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

``` r
library(grid)
library(gtable)
source('biodata/funciones.R')
source('https://raw.githubusercontent.com/maestria-geotel-master/unidad-3-asignacion-1-vecindad-autocorrelacion-espacial/master/lisaclusters.R')
```

### Cargar datos

``` r
load('biodata/Apocynaceae-Meliaceae-Sapotaceae.Rdata')
load('biodata/matriz_ambiental.Rdata')
mi_fam <- mc_apcyn_melic_saptc
mi_fam %>% tibble
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

``` r
bci_env_grid %>% tibble
```

    ## # A tibble: 50 x 39
    ##       id categoria_de_ed… geologia habitat quebrada heterogeneidad_… UTM.EW
    ##    <dbl> <fct>            <fct>    <fct>   <fct>               <dbl>  <dbl>
    ##  1     1 c3               Tb       OldSlo… Yes                0.627  6.26e5
    ##  2     2 c3               Tb       OldLow  Yes                0.394  6.26e5
    ##  3     3 c3               Tb       OldLow  No                 0      6.26e5
    ##  4     4 c3               Tb       OldLow  No                 0      6.26e5
    ##  5     5 c3               Tb       OldSlo… No                 0.461  6.26e5
    ##  6     6 c3               Tb       OldLow  No                 0.0768 6.26e5
    ##  7     7 c3               Tb       OldLow  Yes                0.381  6.26e5
    ##  8     8 c3               Tb       OldLow  Yes                0.211  6.26e5
    ##  9     9 c3               Tb       OldLow  No                 0      6.26e5
    ## 10    10 c3               Tb       OldLow  No                 0      6.26e5
    ## # … with 40 more rows, and 32 more variables: UTM.NS <dbl>,
    ## #   geomorf_llanura_pct <dbl>, geomorf_pico_pct <dbl>,
    ## #   geomorf_interfluvio_pct <dbl>, geomorf_hombrera_pct <dbl>,
    ## #   `geomorf_espolón/gajo_pct` <dbl>, geomorf_vertiente_pct <dbl>,
    ## #   geomorf_vaguada_pct <dbl>, geomorf_piedemonte_pct <dbl>,
    ## #   geomorf_valle_pct <dbl>, geomorf_sima_pct <dbl>, Al <dbl>, B <dbl>,
    ## #   Ca <dbl>, Cu <dbl>, Fe <dbl>, K <dbl>, Mg <dbl>, Mn <dbl>, P <dbl>,
    ## #   Zn <dbl>, N <dbl>, N.min. <dbl>, pH <dbl>, elevacion_media <dbl>,
    ## #   pendiente_media <dbl>, orientacion_media <dbl>,
    ## #   curvatura_perfil_media <dbl>, curvatura_tangencial_media <dbl>,
    ## #   geometry <POLYGON [m]>, abundancia_global <dbl>, riqueza_global <int>

## Preparar datos

### Generar matriz Hellinger

``` r
mi_fam_hel <- decostand (mi_fam, "hellinger")
```

### Transformar matriz ambiental en objeto `sp`, generar vecindad

``` r
bci_env_grid_sp <- bci_env_grid %>% as_Spatial
centroides <- bci_env_grid %>% st_centroid
```

    ## Warning in st_centroid.sf(.): st_centroid assumes attributes are constant
    ## over geometries of x

``` r
bci_xy <- centroides %>% st_coordinates %>% as.data.frame
(vecindad <- bci_env_grid_sp %>% poly2nb)
```

    ## Neighbour list object:
    ## Number of regions: 50 
    ## Number of nonzero links: 314 
    ## Percentage nonzero weights: 12.56 
    ## Average number of links: 6.28

``` r
(pesos_b <- nb2listw(vecindad, style = 'B'))
```

    ## Characteristics of weights list object:
    ## Neighbour list object:
    ## Number of regions: 50 
    ## Number of nonzero links: 314 
    ## Percentage nonzero weights: 12.56 
    ## Average number of links: 6.28 
    ## 
    ## Weights style: B 
    ## Weights constants summary:
    ##    n   nn  S0  S1   S2
    ## B 50 2500 314 628 8488

``` r
plot(bci_env_grid_sp)
plot(vecindad, coords = bci_xy, add=T, col = 'red')
```

![](ee_ecologia_espacial_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

## Autocorrelación espacial

### Autocorrelación espacial de una variable ambiental

``` r
var_ph <- bci_env_grid %>% st_drop_geometry %>% pull(pH)
ph_correl <- sp.correlogram(vecindad,
                            var_ph,
                            order = 9,
                            method = "I",
                            zero.policy = TRUE)
print(ph_correl, p.adj.method = 'holm')
```

    ## Spatial correlogram for var_ph 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.7206808  -0.0204082  0.0060009           9.5667       < 2.2e-16
    ## 2 (50)  0.5123785  -0.0204082  0.0036285           8.8448       < 2.2e-16
    ## 3 (50)  0.2099793  -0.0204082  0.0033715           3.9678       0.0002902
    ## 4 (50) -0.3586281  -0.0204082  0.0044141          -5.0907       1.784e-06
    ## 5 (50) -0.6585951  -0.0204082  0.0073237          -7.4573       6.182e-13
    ## 6 (40) -0.6618167  -0.0256410  0.0089300          -6.7321       1.003e-10
    ## 7 (30) -0.3727208  -0.0344828  0.0113908          -3.1692       0.0045864
    ## 8 (20) -0.2124804  -0.0526316  0.0154419          -1.2864       0.3966415
    ## 9 (10) -0.1056156  -0.1111111  0.0193348           0.0395       0.9684743
    ##           
    ## 1 (50) ***
    ## 2 (50) ***
    ## 3 (50) ***
    ## 4 (50) ***
    ## 5 (50) ***
    ## 6 (40) ***
    ## 7 (30) ** 
    ## 8 (20)    
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(ph_correl)
```

![](ee_ecologia_espacial_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

### Autocorrelación espacial de múltiples variables

#### Autocorrelación espacial de especies (matriz de comunidad)

``` r
suppressWarnings(auto_spp_hel <- calcular_autocorrelacion(
  df_fuente = mi_fam_hel,
  orden = 9,
  obj_vecindad = vecindad,
  pos_var = '(matriz Hellinger)'))
print(auto_spp_hel, p.adj.method = 'holm')
```

    ## $`Aspidosperma spruceanum`
    ## Spatial correlogram for Aspidosperma spruceanum (matriz Hellinger) 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.2821842  -0.0204082  0.0059367           3.9272       0.0007734
    ## 2 (50)  0.0673763  -0.0204082  0.0035903           1.4650       1.0000000
    ## 3 (50) -0.0436623  -0.0204082  0.0033360          -0.4026       1.0000000
    ## 4 (50)  0.0026160  -0.0204082  0.0043679           0.3484       1.0000000
    ## 5 (50) -0.0806656  -0.0204082  0.0072450          -0.7079       1.0000000
    ## 6 (40) -0.1689479  -0.0256410  0.0088074          -1.5270       1.0000000
    ## 7 (30) -0.1566443  -0.0344828  0.0111743          -1.1556       1.0000000
    ## 8 (20) -0.0544376  -0.0526316  0.0149656          -0.0148       1.0000000
    ## 9 (10) -0.0011661  -0.1111111  0.0177650           0.8249       1.0000000
    ##           
    ## 1 (50) ***
    ## 2 (50)    
    ## 3 (50)    
    ## 4 (50)    
    ## 5 (50)    
    ## 6 (40)    
    ## 7 (30)    
    ## 8 (20)    
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`Cedrela odorata`
    ## Spatial correlogram for Cedrela odorata (matriz Hellinger) 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.0680356  -0.0204082  0.0053870           1.2050          1.0000
    ## 2 (50) -0.0263942  -0.0204082  0.0032630          -0.1048          1.0000
    ## 3 (50) -0.0404889  -0.0204082  0.0030321          -0.3647          1.0000
    ## 4 (50) -0.1290646  -0.0204082  0.0039722          -1.7240          0.7624
    ## 5 (50) -0.0233245  -0.0204082  0.0065710          -0.0360          1.0000
    ## 6 (40) -0.0298643  -0.0256410  0.0077574          -0.0480          1.0000
    ## 7 (30)  0.0636205  -0.0344828  0.0093206           1.0162          1.0000
    ## 8 (20)  0.0263009  -0.0526316  0.0108880           0.7565          1.0000
    ## 9 (10) -0.0130695  -0.1111111  0.0043257           1.4907          1.0000
    ## 
    ## $`Chrysophyllum argenteum`
    ## Spatial correlogram for Chrysophyllum argenteum (matriz Hellinger) 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.5664805  -0.0204082  0.0058409           7.6792       1.441e-13
    ## 2 (50)  0.1471605  -0.0204082  0.0035332           2.8191         0.03371
    ## 3 (50) -0.0141407  -0.0204082  0.0032830           0.1094         1.00000
    ## 4 (50) -0.1604742  -0.0204082  0.0042989          -2.1362         0.16330
    ## 5 (50) -0.4333178  -0.0204082  0.0071275          -4.8909       8.031e-06
    ## 6 (40) -0.2670386  -0.0256410  0.0086243          -2.5994         0.05603
    ## 7 (30) -0.0855181  -0.0344828  0.0108512          -0.4899         1.00000
    ## 8 (20) -0.0165061  -0.0526316  0.0142548           0.3026         1.00000
    ## 9 (10) -0.0488776  -0.1111111  0.0154222           0.5011         1.00000
    ##           
    ## 1 (50) ***
    ## 2 (50) *  
    ## 3 (50)    
    ## 4 (50)    
    ## 5 (50) ***
    ## 6 (40) .  
    ## 7 (30)    
    ## 8 (20)    
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`Chrysophyllum cainito`
    ## Spatial correlogram for Chrysophyllum cainito (matriz Hellinger) 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.1708336  -0.0204082  0.0059039           2.4889         0.10250
    ## 2 (50)  0.0955670  -0.0204082  0.0035708           1.9408         0.36596
    ## 3 (50) -0.2002163  -0.0204082  0.0033178          -3.1216         0.01619
    ## 4 (50) -0.0578177  -0.0204082  0.0043443          -0.5676         1.00000
    ## 5 (50) -0.1796886  -0.0204082  0.0072048          -1.8765         0.36596
    ## 6 (40) -0.0394236  -0.0256410  0.0087447          -0.1474         1.00000
    ## 7 (30)  0.0047633  -0.0344828  0.0110636           0.3731         1.00000
    ## 8 (20)  0.0624806  -0.0526316  0.0147221           0.9487         1.00000
    ## 9 (10)  0.0325260  -0.1111111  0.0169624           1.1029         1.00000
    ##         
    ## 1 (50)  
    ## 2 (50)  
    ## 3 (50) *
    ## 4 (50)  
    ## 5 (50)  
    ## 6 (40)  
    ## 7 (30)  
    ## 8 (20)  
    ## 9 (10)  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`Guarea bullata`
    ## Spatial correlogram for Guarea bullata (matriz Hellinger) 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.1510282  -0.0204082  0.0058144           2.2483           0.221
    ## 2 (50)  0.0327022  -0.0204082  0.0035175           0.8955           1.000
    ## 3 (50) -0.0376574  -0.0204082  0.0032684          -0.3017           1.000
    ## 4 (50) -0.0648599  -0.0204082  0.0042799          -0.6795           1.000
    ## 5 (50)  0.0880363  -0.0204082  0.0070950           1.2874           1.000
    ## 6 (40) -0.1562830  -0.0256410  0.0085737          -1.4109           1.000
    ## 7 (30) -0.0344513  -0.0344828  0.0107619           0.0003           1.000
    ## 8 (20) -0.0028253  -0.0526316  0.0140584           0.4201           1.000
    ## 9 (10) -0.1017589  -0.1111111  0.0147748           0.0769           1.000
    ## 
    ## $`Guarea grandifolia`
    ## Spatial correlogram for Guarea grandifolia (matriz Hellinger) 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.0585188  -0.0204082  0.0059923           1.0196               1
    ## 2 (50)  0.0126751  -0.0204082  0.0036234           0.5496               1
    ## 3 (50) -0.0554141  -0.0204082  0.0033667          -0.6033               1
    ## 4 (50)  0.0161066  -0.0204082  0.0044080           0.5500               1
    ## 5 (50) -0.0292339  -0.0204082  0.0073132          -0.1032               1
    ## 6 (40) -0.1143774  -0.0256410  0.0089136          -0.9399               1
    ## 7 (30)  0.0429426  -0.0344828  0.0113619           0.7264               1
    ## 8 (20) -0.0430178  -0.0526316  0.0153782           0.0775               1
    ## 9 (10) -0.0302082  -0.1111111  0.0191248           0.5850               1
    ## 
    ## $`Guarea guidonia`
    ## Spatial correlogram for Guarea guidonia (matriz Hellinger) 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.2292987  -0.0204082  0.0059067           3.2491         0.01042
    ## 2 (50)  0.0958004  -0.0204082  0.0035724           1.9443         0.36304
    ## 3 (50) -0.1561422  -0.0204082  0.0033194          -2.3559         0.14782
    ## 4 (50) -0.1444407  -0.0204082  0.0043463          -1.8814         0.36304
    ## 5 (50) -0.0247025  -0.0204082  0.0072082          -0.0506         1.00000
    ## 6 (40)  0.0659767  -0.0256410  0.0087501           0.9794         1.00000
    ## 7 (30)  0.0277528  -0.0344828  0.0110731           0.5914         1.00000
    ## 8 (20) -0.1975110  -0.0526316  0.0147430          -1.1932         1.00000
    ## 9 (10)  0.0552299  -0.1111111  0.0170315           1.2746         1.00000
    ##         
    ## 1 (50) *
    ## 2 (50)  
    ## 3 (50)  
    ## 4 (50)  
    ## 5 (50)  
    ## 6 (40)  
    ## 7 (30)  
    ## 8 (20)  
    ## 9 (10)  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`Lacmellea panamensis`
    ## Spatial correlogram for Lacmellea panamensis (matriz Hellinger) 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.1321626  -0.0204082  0.0059197           1.9830          0.4263
    ## 2 (50)  0.0473699  -0.0204082  0.0035802           1.1328          1.0000
    ## 3 (50) -0.0128834  -0.0204082  0.0033266           0.1305          1.0000
    ## 4 (50) -0.0364279  -0.0204082  0.0043557          -0.2427          1.0000
    ## 5 (50) -0.1073805  -0.0204082  0.0072242          -1.0233          1.0000
    ## 6 (40) -0.1054048  -0.0256410  0.0087750          -0.8515          1.0000
    ## 7 (30) -0.0907490  -0.0344828  0.0111171          -0.5336          1.0000
    ## 8 (20) -0.0528651  -0.0526316  0.0148398          -0.0019          1.0000
    ## 9 (10) -0.0067893  -0.1111111  0.0173505           0.7920          1.0000
    ## 
    ## $`Pouteria fossicola`
    ## Spatial correlogram for Pouteria fossicola (matriz Hellinger) 
    ## method: Moran's I

    ## Warning in sqrt(res[, 3]): NaNs produced

    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50) -0.0529694  -0.0204082  0.0042919          -0.4970               1
    ## 2 (50) -0.0453485  -0.0204082  0.0026109          -0.4881               1
    ## 3 (50) -0.0327410  -0.0204082  0.0024268          -0.2504               1
    ## 4 (50)  0.0462460  -0.0204082  0.0031839           1.1813               1
    ## 5 (50) -0.0633712  -0.0204082  0.0052283          -0.5942               1
    ## 6 (40)  0.0044916  -0.0256410  0.0056657           0.4003               1
    ## 7 (30) -0.0081827  -0.0344828  0.0056278           0.3506               1
    ## 8 (20) -0.0208569  -0.0526316  0.0027650           0.6043               1
    ## 9 (10) -0.0335312  -0.1111111 -0.0224473               NA              NA
    ## 
    ## $`Pouteria reticulata`
    ## Spatial correlogram for Pouteria reticulata (matriz Hellinger) 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.2969399  -0.0204082  0.0057984           4.1676       0.0002771
    ## 2 (50)  0.2170940  -0.0204082  0.0035079           4.0100       0.0004858
    ## 3 (50)  0.0876786  -0.0204082  0.0032595           1.8932       0.2916631
    ## 4 (50) -0.1432722  -0.0204082  0.0042684          -1.8806       0.2916631
    ## 5 (50) -0.2498965  -0.0204082  0.0070754          -2.7282       0.0445700
    ## 6 (40) -0.2386480  -0.0256410  0.0085432          -2.3045       0.1271549
    ## 7 (30) -0.1843909  -0.0344828  0.0107079          -1.4487       0.4422814
    ## 8 (20) -0.1238287  -0.0526316  0.0139397          -0.6030       0.9117406
    ## 9 (10) -0.0216820  -0.1111111  0.0143837           0.7457       0.9117406
    ##           
    ## 1 (50) ***
    ## 2 (50) ***
    ## 3 (50)    
    ## 4 (50)    
    ## 5 (50) *  
    ## 6 (40)    
    ## 7 (30)    
    ## 8 (20)    
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`Pouteria stipitata`
    ## Spatial correlogram for Pouteria stipitata (matriz Hellinger) 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.1227159  -0.0204082  0.0058932           1.8644         0.37361
    ## 2 (50)  0.0989802  -0.0204082  0.0035644           1.9997         0.31872
    ## 3 (50) -0.0105557  -0.0204082  0.0033120           0.1712         1.00000
    ## 4 (50) -0.1673297  -0.0204082  0.0043366          -2.2310         0.20542
    ## 5 (50) -0.2540736  -0.0204082  0.0071917          -2.7554         0.05277
    ## 6 (40) -0.1161225  -0.0256410  0.0087244          -0.9687         1.00000
    ## 7 (30) -0.0952540  -0.0344828  0.0110278          -0.5787         1.00000
    ## 8 (20)  0.0584658  -0.0526316  0.0146433           0.9181         1.00000
    ## 9 (10)  0.0094232  -0.1111111  0.0167026           0.9326         1.00000
    ##         
    ## 1 (50)  
    ## 2 (50)  
    ## 3 (50)  
    ## 4 (50)  
    ## 5 (50) .
    ## 6 (40)  
    ## 7 (30)  
    ## 8 (20)  
    ## 9 (10)  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`Rauvolfia littoralis`
    ## Spatial correlogram for Rauvolfia littoralis (matriz Hellinger) 
    ## method: Moran's I

    ## Warning in sqrt(res[, 3]): NaNs produced

    ##           estimate expectation    variance standard deviate
    ## 1 (50) -3.0782e-02 -2.0408e-02  2.7954e-05          -1.9622
    ## 2 (50) -1.5426e-02 -2.0408e-02  7.1776e-05           0.5881
    ## 3 (50) -1.7206e-02 -2.0408e-02  6.9697e-05           0.3836
    ## 4 (50) -8.5234e-03 -2.0408e-02  1.1427e-04           1.1118
    ## 5 (50) -2.0408e-02 -2.0408e-02  4.3368e-19           0.0000
    ## 6 (40) -2.4490e-02 -2.5641e-02 -2.4789e-03               NA
    ## 7 (30) -2.8571e-02 -3.4483e-02 -8.7514e-03               NA
    ## 8 (20) -3.2653e-02 -5.2632e-02 -2.8865e-02               NA
    ## 9 (10)  4.0816e-03 -1.1111e-01 -1.2670e-01               NA
    ##        Pr(I) two sided
    ## 1 (50)          0.2487
    ## 2 (50)          1.0000
    ## 3 (50)          1.0000
    ## 4 (50)          1.0000
    ## 5 (50)          1.0000
    ## 6 (40)              NA
    ## 7 (30)              NA
    ## 8 (20)              NA
    ## 9 (10)              NA
    ## 
    ## $`Tabernaemontana arborea`
    ## Spatial correlogram for Tabernaemontana arborea (matriz Hellinger) 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.6834557  -0.0204082  0.0060406           9.0563       < 2.2e-16
    ## 2 (50)  0.5498674  -0.0204082  0.0036522           9.4365       < 2.2e-16
    ## 3 (50)  0.3583748  -0.0204082  0.0033934           6.5024       3.953e-10
    ## 4 (50) -0.1243569  -0.0204082  0.0044427          -1.5595         0.23774
    ## 5 (50) -0.6488906  -0.0204082  0.0073724          -7.3196       1.741e-12
    ## 6 (40) -0.6762620  -0.0256410  0.0090058          -6.8559       4.251e-11
    ## 7 (30) -0.5063934  -0.0344828  0.0115247          -4.3959       4.413e-05
    ## 8 (20) -0.3843947  -0.0526316  0.0157363          -2.6447         0.02453
    ## 9 (10) -0.2253977  -0.1111111  0.0203053          -0.8020         0.42254
    ##           
    ## 1 (50) ***
    ## 2 (50) ***
    ## 3 (50) ***
    ## 4 (50)    
    ## 5 (50) ***
    ## 6 (40) ***
    ## 7 (30) ***
    ## 8 (20) *  
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`Thevetia ahouai`
    ## Spatial correlogram for Thevetia ahouai (matriz Hellinger) 
    ## method: Moran's I

    ## Warning in sqrt(res[, 3]): NaNs produced

    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.0073827  -0.0204082  0.0050239           0.3921          1.0000
    ## 2 (50) -0.1168566  -0.0204082  0.0030467          -1.7473          0.6446
    ## 3 (50)  0.0020077  -0.0204082  0.0028314           0.4213          1.0000
    ## 4 (50) -0.0273975  -0.0204082  0.0037108          -0.1147          1.0000
    ## 5 (50)  0.0515544  -0.0204082  0.0061257           0.9194          1.0000
    ## 6 (40) -0.0238811  -0.0256410  0.0070638           0.0209          1.0000
    ## 7 (30) -0.0842869  -0.0344828  0.0080960          -0.5535          1.0000
    ## 8 (20)  0.0689993  -0.0526316  0.0081944           1.3437          1.0000
    ## 9 (10)  0.0031630  -0.1111111 -0.0045524               NA              NA
    ## 
    ## $`Trichilia pallida`
    ## Spatial correlogram for Trichilia pallida (matriz Hellinger) 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.4302533  -0.0204082  0.0059272           5.8537       4.328e-08
    ## 2 (50)  0.2293836  -0.0204082  0.0035846           4.1721       0.0002112
    ## 3 (50)  0.0225725  -0.0204082  0.0033307           0.7447       1.0000000
    ## 4 (50)  0.0058207  -0.0204082  0.0043611           0.3972       1.0000000
    ## 5 (50) -0.1947077  -0.0204082  0.0072333          -2.0494       0.2021126
    ## 6 (40) -0.4666369  -0.0256410  0.0087891          -4.7039       2.042e-05
    ## 7 (30) -0.3298341  -0.0344828  0.0111421          -2.7980       0.0308476
    ## 8 (20) -0.0940116  -0.0526316  0.0148948          -0.3391       1.0000000
    ## 9 (10)  0.0031417  -0.1111111  0.0175316           0.8629       1.0000000
    ##           
    ## 1 (50) ***
    ## 2 (50) ***
    ## 3 (50)    
    ## 4 (50)    
    ## 5 (50)    
    ## 6 (40) ***
    ## 7 (30) *  
    ## 8 (20)    
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`Trichilia tuberculata`
    ## Spatial correlogram for Trichilia tuberculata (matriz Hellinger) 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.5819564  -0.0204082  0.0059458           7.8118       5.072e-14
    ## 2 (50)  0.4050703  -0.0204082  0.0035957           7.0955       1.031e-11
    ## 3 (50)  0.0995892  -0.0204082  0.0033410           2.0760          0.1137
    ## 4 (50) -0.1726643  -0.0204082  0.0043745          -2.3020          0.1067
    ## 5 (50) -0.4885437  -0.0204082  0.0072562          -5.4956       2.725e-07
    ## 6 (40) -0.4566198  -0.0256410  0.0088248          -4.5878       2.688e-05
    ## 7 (30) -0.2672493  -0.0344828  0.0112051          -2.1989          0.1115
    ## 8 (20) -0.1658190  -0.0526316  0.0150332          -0.9231          0.7119
    ## 9 (10) -0.0928353  -0.1111111  0.0179879           0.1363          0.8916
    ##           
    ## 1 (50) ***
    ## 2 (50) ***
    ## 3 (50)    
    ## 4 (50)    
    ## 5 (50) ***
    ## 6 (40) ***
    ## 7 (30)    
    ## 8 (20)    
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
dim_panel <- rev(n2mfrow(ncol(mi_fam_hel)))
```

``` r
par(mfrow = dim_panel)
suppressWarnings(invisible(lapply(auto_spp_hel, function(x) plot(x, main = x$var))))
```

![](ee_ecologia_espacial_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

#### Autocorrelación espacial de datos ambientales (matriz ambiental)

``` r
bci_env_grid_num <- bci_env_grid %>%
  st_drop_geometry %>% 
  select_if(is.numeric) %>% 
  select(-id, -UTM.EW, -UTM.NS)
suppressWarnings(auto_amb <- calcular_autocorrelacion(
  df_fuente = bci_env_grid_num,
  orden = 9,
  obj_vecindad = vecindad))
print(auto_amb, p.adj.method = 'holm')
```

    ## $heterogeneidad_ambiental
    ## Spatial correlogram for heterogeneidad_ambiental 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.0860024  -0.0204082  0.0060473           1.3684               1
    ## 2 (50) -0.0627458  -0.0204082  0.0036561          -0.7002               1
    ## 3 (50) -0.0645692  -0.0204082  0.0033971          -0.7577               1
    ## 4 (50) -0.0084073  -0.0204082  0.0044475           0.1800               1
    ## 5 (50)  0.0423991  -0.0204082  0.0073806           0.7311               1
    ## 6 (40)  0.0099421  -0.0256410  0.0090185           0.3747               1
    ## 7 (30) -0.0221393  -0.0344828  0.0115472           0.1149               1
    ## 8 (20) -0.0706705  -0.0526316  0.0157857          -0.1436               1
    ## 9 (10) -0.0057738  -0.1111111  0.0204681           0.7363               1
    ## 
    ## $geomorf_llanura_pct
    ## Spatial correlogram for geomorf_llanura_pct 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.4143623  -0.0204082  0.0055958           5.8121       5.553e-08
    ## 2 (50) -0.0367893  -0.0204082  0.0033873          -0.2815          1.0000
    ## 3 (50) -0.1372856  -0.0204082  0.0031475          -2.0833          0.2978
    ## 4 (50) -0.0849775  -0.0204082  0.0041225          -1.0056          1.0000
    ## 5 (50) -0.0473605  -0.0204082  0.0068270          -0.3262          1.0000
    ## 6 (40) -0.1158353  -0.0256410  0.0081561          -0.9987          1.0000
    ## 7 (30)  0.0014408  -0.0344828  0.0100246           0.3588          1.0000
    ## 8 (20)  0.0643630  -0.0526316  0.0124365           1.0491          1.0000
    ## 9 (10)  0.0175592  -0.1111111  0.0094293           1.3251          1.0000
    ##           
    ## 1 (50) ***
    ## 2 (50)    
    ## 3 (50)    
    ## 4 (50)    
    ## 5 (50)    
    ## 6 (40)    
    ## 7 (30)    
    ## 8 (20)    
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $geomorf_pico_pct
    ## Spatial correlogram for geomorf_pico_pct 
    ## method: Moran's I

    ## Warning in sqrt(res[, 3]): NaNs produced

    ##           estimate expectation    variance standard deviate
    ## 1 (50)  0.03862269 -0.02040816  0.00141273           1.5705
    ## 2 (50) -0.03069198 -0.02040816  0.00089638          -0.3435
    ## 3 (50) -0.03853384 -0.02040816  0.00083519          -0.6272
    ## 4 (50)  0.03523986 -0.02040816  0.00111116           1.6694
    ## 5 (50) -0.04169381 -0.02040816  0.00169795          -0.5166
    ## 6 (40) -0.03960912 -0.02564103  0.00016615          -1.0837
    ## 7 (30) -0.04794788 -0.03448276 -0.00408156               NA
    ## 8 (20) -0.05628664 -0.05263158 -0.01859262               NA
    ## 9 (10)  0.00833876 -0.11111111 -0.09284012               NA
    ##        Pr(I) two sided
    ## 1 (50)          0.5814
    ## 2 (50)          1.0000
    ## 3 (50)          1.0000
    ## 4 (50)          0.5702
    ## 5 (50)          1.0000
    ## 6 (40)          1.0000
    ## 7 (30)              NA
    ## 8 (20)              NA
    ## 9 (10)              NA
    ## 
    ## $geomorf_interfluvio_pct
    ## Spatial correlogram for geomorf_interfluvio_pct 
    ## method: Moran's I

    ## Warning in sqrt(res[, 3]): NaNs produced

    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.1103328  -0.0204082  0.0047651           1.8940          0.4658
    ## 2 (50) -0.0113988  -0.0204082  0.0028926           0.1675          1.0000
    ## 3 (50) -0.1040701  -0.0204082  0.0026883          -1.6136          0.7464
    ## 4 (50) -0.0361553  -0.0204082  0.0035245          -0.2652          1.0000
    ## 5 (50) -0.0131219  -0.0204082  0.0058084           0.0956          1.0000
    ## 6 (40)  0.0176134  -0.0256410  0.0065695           0.5337          1.0000
    ## 7 (30) -0.0026818  -0.0344828  0.0072234           0.3742          1.0000
    ## 8 (20) -0.0027238  -0.0526316  0.0062747           0.6300          1.0000
    ## 9 (10)  0.0046557  -0.1111111 -0.0108795               NA              NA
    ## 
    ## $geomorf_hombrera_pct
    ## Spatial correlogram for geomorf_hombrera_pct 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.1432057  -0.0204082  0.0058532           2.1386          0.2922
    ## 2 (50)  0.0337750  -0.0204082  0.0035406           0.9106          1.0000
    ## 3 (50) -0.0617868  -0.0204082  0.0032898          -0.7214          1.0000
    ## 4 (50) -0.0404654  -0.0204082  0.0043078          -0.3056          1.0000
    ## 5 (50)  0.0073955  -0.0204082  0.0071426           0.3290          1.0000
    ## 6 (40) -0.0262658  -0.0256410  0.0086479          -0.0067          1.0000
    ## 7 (30) -0.0287347  -0.0344828  0.0108928           0.0551          1.0000
    ## 8 (20) -0.0073196  -0.0526316  0.0143463           0.3783          1.0000
    ## 9 (10) -0.1107890  -0.1111111  0.0157237           0.0026          1.0000
    ## 
    ## $`geomorf_espolón/gajo_pct`
    ## Spatial correlogram for geomorf_espolón/gajo_pct 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.3144587  -0.0204082  0.0058301           4.3857       0.0001041
    ## 2 (50)  0.0621159  -0.0204082  0.0035268           1.3896       1.0000000
    ## 3 (50) -0.0852734  -0.0204082  0.0032771          -1.1331       1.0000000
    ## 4 (50) -0.1559902  -0.0204082  0.0042912          -2.0697       0.3078214
    ## 5 (50) -0.0631538  -0.0204082  0.0071143          -0.5068       1.0000000
    ## 6 (40)  0.0377016  -0.0256410  0.0086037           0.6829       1.0000000
    ## 7 (30)  0.0377135  -0.0344828  0.0108148           0.6942       1.0000000
    ## 8 (20) -0.0486254  -0.0526316  0.0141748           0.0336       1.0000000
    ## 9 (10) -0.0120589  -0.1111111  0.0151587           0.8045       1.0000000
    ##           
    ## 1 (50) ***
    ## 2 (50)    
    ## 3 (50)    
    ## 4 (50)    
    ## 5 (50)    
    ## 6 (40)    
    ## 7 (30)    
    ## 8 (20)    
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $geomorf_vertiente_pct
    ## Spatial correlogram for geomorf_vertiente_pct 
    ## method: Moran's I
    ##           estimate expectation    variance standard deviate
    ## 1 (50)  0.29100690 -0.02040816  0.00571074           4.1209
    ## 2 (50) -0.08242984 -0.02040816  0.00345575          -1.0550
    ## 3 (50) -0.04147761 -0.02040816  0.00321109          -0.3718
    ## 4 (50) -0.03722678 -0.02040816  0.00420526          -0.2594
    ## 5 (50) -0.11889259 -0.02040816  0.00696796          -1.1798
    ## 6 (40) -0.12709127 -0.02564103  0.00837576          -1.1085
    ## 7 (30)  0.03934544 -0.03448276  0.01041234           0.7235
    ## 8 (20)  0.03963128 -0.05263158  0.01328947           0.8003
    ## 9 (10)  0.00097986 -0.11111111  0.01224063           1.0131
    ##        Pr(I) two sided    
    ## 1 (50)       0.0003396 ***
    ## 2 (50)       1.0000000    
    ## 3 (50)       1.0000000    
    ## 4 (50)       1.0000000    
    ## 5 (50)       1.0000000    
    ## 6 (40)       1.0000000    
    ## 7 (30)       1.0000000    
    ## 8 (20)       1.0000000    
    ## 9 (10)       1.0000000    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $geomorf_vaguada_pct
    ## Spatial correlogram for geomorf_vaguada_pct 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.3820769  -0.0204082  0.0057277           5.3181       9.436e-07
    ## 2 (50)  0.0835405  -0.0204082  0.0034659           1.7657         0.46469
    ## 3 (50) -0.1747488  -0.0204082  0.0032205          -2.7197         0.05227
    ## 4 (50) -0.1680424  -0.0204082  0.0042175          -2.2733         0.16105
    ## 5 (50) -0.0172009  -0.0204082  0.0069888           0.0384         1.00000
    ## 6 (40)  0.0039059  -0.0256410  0.0084082           0.3222         1.00000
    ## 7 (30) -0.0150833  -0.0344828  0.0104696           0.1896         1.00000
    ## 8 (20) -0.0399341  -0.0526316  0.0134155           0.1096         1.00000
    ## 9 (10) -0.0179433  -0.1111111  0.0126559           0.8282         1.00000
    ##           
    ## 1 (50) ***
    ## 2 (50)    
    ## 3 (50) .  
    ## 4 (50)    
    ## 5 (50)    
    ## 6 (40)    
    ## 7 (30)    
    ## 8 (20)    
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $geomorf_piedemonte_pct
    ## Spatial correlogram for geomorf_piedemonte_pct 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50) -0.0380217  -0.0204082  0.0054299          -0.2390          1.0000
    ## 2 (50) -0.0433465  -0.0204082  0.0032885          -0.4000          1.0000
    ## 3 (50)  0.1007212  -0.0204082  0.0030559           2.1912          0.2559
    ## 4 (50) -0.1523801  -0.0204082  0.0040031          -2.0858          0.2959
    ## 5 (50) -0.0515249  -0.0204082  0.0066237          -0.3823          1.0000
    ## 6 (40)  0.0169652  -0.0256410  0.0078394           0.4812          1.0000
    ## 7 (30)  0.0189232  -0.0344828  0.0094654           0.5489          1.0000
    ## 8 (20)  0.0211800  -0.0526316  0.0112065           0.6973          1.0000
    ## 9 (10)  0.0069996  -0.1111111  0.0053754           1.6110          0.7503
    ## 
    ## $geomorf_valle_pct
    ## Spatial correlogram for geomorf_valle_pct 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.1019140  -0.0204082  0.0054158           1.6622          0.7719
    ## 2 (50)  0.0425053  -0.0204082  0.0032801           1.0985          1.0000
    ## 3 (50)  0.0139414  -0.0204082  0.0030481           0.6222          1.0000
    ## 4 (50) -0.0635342  -0.0204082  0.0039930          -0.6825          1.0000
    ## 5 (50) -0.0758324  -0.0204082  0.0066064          -0.6819          1.0000
    ## 6 (40) -0.0828752  -0.0256410  0.0078125          -0.6475          1.0000
    ## 7 (30) -0.1017703  -0.0344828  0.0094179          -0.6934          1.0000
    ## 8 (20) -0.0430101  -0.0526316  0.0111020           0.0913          1.0000
    ## 9 (10)  0.0118862  -0.1111111  0.0050307           1.7341          0.7461
    ## 
    ## $geomorf_sima_pct
    ## Spatial correlogram for geomorf_sima_pct 
    ## method: Moran's I

    ## Warning in sqrt(res[, 3]): NaNs produced

    ##           estimate expectation    variance standard deviate
    ## 1 (50) -3.0782e-02 -2.0408e-02  2.7954e-05          -1.9622
    ## 2 (50) -1.5426e-02 -2.0408e-02  7.1776e-05           0.5881
    ## 3 (50) -1.7206e-02 -2.0408e-02  6.9697e-05           0.3836
    ## 4 (50) -8.5234e-03 -2.0408e-02  1.1427e-04           1.1118
    ## 5 (50) -2.0408e-02 -2.0408e-02  3.7947e-18           0.0000
    ## 6 (40) -2.4490e-02 -2.5641e-02 -2.4789e-03               NA
    ## 7 (30) -2.8571e-02 -3.4483e-02 -8.7514e-03               NA
    ## 8 (20) -3.2653e-02 -5.2632e-02 -2.8865e-02               NA
    ## 9 (10)  4.0816e-03 -1.1111e-01 -1.2670e-01               NA
    ##        Pr(I) two sided
    ## 1 (50)          0.2487
    ## 2 (50)          1.0000
    ## 3 (50)          1.0000
    ## 4 (50)          1.0000
    ## 5 (50)          1.0000
    ## 6 (40)              NA
    ## 7 (30)              NA
    ## 8 (20)              NA
    ## 9 (10)              NA
    ## 
    ## $Al
    ## Spatial correlogram for Al 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.3818824  -0.0204082  0.0059634           5.2095       1.704e-06
    ## 2 (50)  0.0570533  -0.0204082  0.0036062           1.2899         0.79574
    ## 3 (50) -0.0709217  -0.0204082  0.0033507          -0.8726         1.00000
    ## 4 (50)  0.0962211  -0.0204082  0.0043871           1.7608         0.46960
    ## 5 (50)  0.0052579  -0.0204082  0.0072777           0.3009         1.00000
    ## 6 (40) -0.2391856  -0.0256410  0.0088583          -2.2689         0.16292
    ## 7 (30) -0.3268705  -0.0344828  0.0112642          -2.7549         0.04697
    ## 8 (20) -0.2260053  -0.0526316  0.0151633          -1.4079         0.79574
    ## 9 (10)  0.0200170  -0.1111111  0.0184167           0.9663         1.00000
    ##           
    ## 1 (50) ***
    ## 2 (50)    
    ## 3 (50)    
    ## 4 (50)    
    ## 5 (50)    
    ## 6 (40)    
    ## 7 (30) *  
    ## 8 (20)    
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $B
    ## Spatial correlogram for B 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.6609907  -0.0204082  0.0059727           8.8169       < 2.2e-16
    ## 2 (50)  0.4605069  -0.0204082  0.0036117           8.0022       9.774e-15
    ## 3 (50)  0.2413760  -0.0204082  0.0033559           4.5190       2.486e-05
    ## 4 (50) -0.0838790  -0.0204082  0.0043938          -0.9575          0.6766
    ## 5 (50) -0.4166034  -0.0204082  0.0072891          -4.6406       1.737e-05
    ## 6 (40) -0.5790424  -0.0256410  0.0088761          -5.8739       2.979e-08
    ## 7 (30) -0.5395948  -0.0344828  0.0112956          -4.7526       1.205e-05
    ## 8 (20) -0.4105206  -0.0526316  0.0152324          -2.8998          0.0112
    ## 9 (10) -0.1525029  -0.1111111  0.0186443          -0.3031          0.7618
    ##           
    ## 1 (50) ***
    ## 2 (50) ***
    ## 3 (50) ***
    ## 4 (50)    
    ## 5 (50) ***
    ## 6 (40) ***
    ## 7 (30) ***
    ## 8 (20) *  
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $Ca
    ## Spatial correlogram for Ca 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.6912755  -0.0204082  0.0058563           9.2998       < 2.2e-16
    ## 2 (50)  0.4427179  -0.0204082  0.0035424           7.7812       5.745e-14
    ## 3 (50)  0.1827536  -0.0204082  0.0032915           3.5411        0.001992
    ## 4 (50) -0.0629086  -0.0204082  0.0043100          -0.6474        1.000000
    ## 5 (50) -0.2654968  -0.0204082  0.0071464          -2.8992        0.011223
    ## 6 (40) -0.5340380  -0.0256410  0.0086538          -5.4651       3.238e-07
    ## 7 (30) -0.5108024  -0.0344828  0.0109032          -4.5617       3.045e-05
    ## 8 (20) -0.4509107  -0.0526316  0.0143692          -3.3226        0.003568
    ## 9 (10) -0.1705730  -0.1111111  0.0157992          -0.4731        1.000000
    ##           
    ## 1 (50) ***
    ## 2 (50) ***
    ## 3 (50) ** 
    ## 4 (50)    
    ## 5 (50) *  
    ## 6 (40) ***
    ## 7 (30) ***
    ## 8 (20) ** 
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $Cu
    ## Spatial correlogram for Cu 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.5226695  -0.0204082  0.0058150           7.1218        9.59e-12
    ## 2 (50)  0.2416759  -0.0204082  0.0035178           4.4188        7.94e-05
    ## 3 (50) -0.0203270  -0.0204082  0.0032687           0.0014         1.00000
    ## 4 (50) -0.1854098  -0.0204082  0.0042803          -2.5220         0.05834
    ## 5 (50) -0.2395438  -0.0204082  0.0070958          -2.6014         0.05570
    ## 6 (40) -0.2745080  -0.0256410  0.0085749          -2.6875         0.05039
    ## 7 (30) -0.2243511  -0.0344828  0.0107639          -1.8301         0.26896
    ## 8 (20) -0.1906877  -0.0526316  0.0140627          -1.1642         0.73305
    ## 9 (10) -0.0720562  -0.1111111  0.0147892           0.3211         1.00000
    ##           
    ## 1 (50) ***
    ## 2 (50) ***
    ## 3 (50)    
    ## 4 (50) .  
    ## 5 (50) .  
    ## 6 (40) .  
    ## 7 (30)    
    ## 8 (20)    
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $Fe
    ## Spatial correlogram for Fe 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.6149322  -0.0204082  0.0058727           8.2906       1.014e-15
    ## 2 (50)  0.2291840  -0.0204082  0.0035522           4.1878       0.0002254
    ## 3 (50) -0.1215160  -0.0204082  0.0033006          -1.7599       0.3137023
    ## 4 (50) -0.0068754  -0.0204082  0.0043219           0.2059       1.0000000
    ## 5 (50)  0.1537142  -0.0204082  0.0071666           2.0568       0.1985109
    ## 6 (40) -0.0779006  -0.0256410  0.0086851          -0.5608       1.0000000
    ## 7 (30) -0.3642717  -0.0344828  0.0109585          -3.1504       0.0097843
    ## 8 (20) -0.4964697  -0.0526316  0.0144910          -3.6870       0.0015882
    ## 9 (10) -0.2123717  -0.1111111  0.0162006          -0.7956       1.0000000
    ##           
    ## 1 (50) ***
    ## 2 (50) ***
    ## 3 (50)    
    ## 4 (50)    
    ## 5 (50)    
    ## 6 (40)    
    ## 7 (30) ** 
    ## 8 (20) ** 
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $K
    ## Spatial correlogram for K 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.7423004  -0.0204082  0.0059452           9.8918       < 2.2e-16
    ## 2 (50)  0.5183397  -0.0204082  0.0035954           8.9849       < 2.2e-16
    ## 3 (50)  0.2195267  -0.0204082  0.0033407           4.1512       0.0001654
    ## 4 (50) -0.0540866  -0.0204082  0.0043740          -0.5092       0.6105940
    ## 5 (50) -0.2426115  -0.0204082  0.0072554          -2.6087       0.0272686
    ## 6 (40) -0.5610304  -0.0256410  0.0088236          -5.6996       8.405e-08
    ## 7 (30) -0.6277187  -0.0344828  0.0112030          -5.6048       1.251e-07
    ## 8 (20) -0.5285017  -0.0526316  0.0150287          -3.8818       0.0004148
    ## 9 (10) -0.2704600  -0.1111111  0.0179729          -1.1886       0.4691840
    ##           
    ## 1 (50) ***
    ## 2 (50) ***
    ## 3 (50) ***
    ## 4 (50)    
    ## 5 (50) *  
    ## 6 (40) ***
    ## 7 (30) ***
    ## 8 (20) ***
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $Mg
    ## Spatial correlogram for Mg 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.6280982  -0.0204082  0.0056852           8.6009       < 2.2e-16
    ## 2 (50)  0.3284847  -0.0204082  0.0034405           5.9481       2.170e-08
    ## 3 (50)  0.1441299  -0.0204082  0.0031970           2.9100       0.0144557
    ## 4 (50) -0.0036289  -0.0204082  0.0041869           0.2593       0.8641870
    ## 5 (50) -0.1246833  -0.0204082  0.0069366          -1.2520       0.6317019
    ## 6 (40) -0.4156643  -0.0256410  0.0083270          -4.2741       0.0001151
    ## 7 (30) -0.4765721  -0.0344828  0.0103262          -4.3505       9.507e-05
    ## 8 (20) -0.4149883  -0.0526316  0.0130999          -3.1659       0.0077294
    ## 9 (10) -0.1957826  -0.1111111  0.0116160          -0.7856       0.8641870
    ##           
    ## 1 (50) ***
    ## 2 (50) ***
    ## 3 (50) *  
    ## 4 (50)    
    ## 5 (50)    
    ## 6 (40) ***
    ## 7 (30) ***
    ## 8 (20) ** 
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $Mn
    ## Spatial correlogram for Mn 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.4428375  -0.0204082  0.0059686           5.9962       1.818e-08
    ## 2 (50)  0.0998225  -0.0204082  0.0036093           2.0013       0.3175487
    ## 3 (50) -0.1231068  -0.0204082  0.0033536          -1.7734       0.4569707
    ## 4 (50) -0.2926423  -0.0204082  0.0043909          -4.1083       0.0003188
    ## 5 (50) -0.1490636  -0.0204082  0.0072841          -1.5074       0.6584911
    ## 6 (40) -0.0306688  -0.0256410  0.0088683          -0.0534       1.0000000
    ## 7 (30) -0.0694226  -0.0344828  0.0112819          -0.3290       1.0000000
    ## 8 (20) -0.0795381  -0.0526316  0.0152022          -0.2182       1.0000000
    ## 9 (10) -0.0235811  -0.1111111  0.0185448           0.6428       1.0000000
    ##           
    ## 1 (50) ***
    ## 2 (50)    
    ## 3 (50)    
    ## 4 (50) ***
    ## 5 (50)    
    ## 6 (40)    
    ## 7 (30)    
    ## 8 (20)    
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $P
    ## Spatial correlogram for P 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.3691771  -0.0204082  0.0058850           5.0784       3.425e-06
    ## 2 (50) -0.1858535  -0.0204082  0.0035595          -2.7731        0.038872
    ## 3 (50) -0.0916332  -0.0204082  0.0033074          -1.2385        0.431077
    ## 4 (50)  0.1552794  -0.0204082  0.0043307           2.6697        0.045552
    ## 5 (50) -0.1475605  -0.0204082  0.0071816          -1.5004        0.400514
    ## 6 (40) -0.3432894  -0.0256410  0.0087086          -3.4039        0.005315
    ## 7 (30) -0.1306802  -0.0344828  0.0109999          -0.9172        0.431077
    ## 8 (20)  0.1819688  -0.0526316  0.0145820           1.9428        0.260221
    ## 9 (10)  0.1326453  -0.1111111  0.0165006           1.8976        0.260221
    ##           
    ## 1 (50) ***
    ## 2 (50) *  
    ## 3 (50)    
    ## 4 (50) *  
    ## 5 (50)    
    ## 6 (40) ** 
    ## 7 (30)    
    ## 8 (20)    
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $Zn
    ## Spatial correlogram for Zn 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.8594295  -0.0204082  0.0059101          11.4447       < 2.2e-16
    ## 2 (50)  0.6434897  -0.0204082  0.0035745          11.1044       < 2.2e-16
    ## 3 (50)  0.3581624  -0.0204082  0.0033213           6.5689       3.041e-10
    ## 4 (50) -0.0232544  -0.0204082  0.0043488          -0.0432          0.9656
    ## 5 (50) -0.4952279  -0.0204082  0.0072124          -5.5910       9.031e-08
    ## 6 (40) -0.6619180  -0.0256410  0.0087565          -6.7996       7.346e-11
    ## 7 (30) -0.6783855  -0.0344828  0.0110846          -6.1159       4.800e-09
    ## 8 (20) -0.6365820  -0.0526316  0.0147681          -4.8052       4.638e-06
    ## 9 (10) -0.3666944  -0.1111111  0.0171142          -1.9537          0.1015
    ##           
    ## 1 (50) ***
    ## 2 (50) ***
    ## 3 (50) ***
    ## 4 (50)    
    ## 5 (50) ***
    ## 6 (40) ***
    ## 7 (30) ***
    ## 8 (20) ***
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $N
    ## Spatial correlogram for N 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.6282989  -0.0204082  0.0060196           8.3611       4.969e-16
    ## 2 (50)  0.3128380  -0.0204082  0.0036397           5.5238       1.991e-07
    ## 3 (50)  0.0973180  -0.0204082  0.0033818           2.0244        0.171713
    ## 4 (50) -0.4825354  -0.0204082  0.0044276          -6.9451       2.648e-11
    ## 5 (50) -0.8074332  -0.0204082  0.0073466          -9.1821       < 2.2e-16
    ## 6 (40) -0.3514819  -0.0256410  0.0089657          -3.4412        0.002895
    ## 7 (30) -0.1627962  -0.0344828  0.0114538          -1.1989        0.461105
    ## 8 (20)  0.0154467  -0.0526316  0.0155805           0.5454        0.585475
    ## 9 (10)  0.1333637  -0.1111111  0.0197916           1.7378        0.246750
    ##           
    ## 1 (50) ***
    ## 2 (50) ***
    ## 3 (50)    
    ## 4 (50) ***
    ## 5 (50) ***
    ## 6 (40) ** 
    ## 7 (30)    
    ## 8 (20)    
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $N.min.
    ## Spatial correlogram for N.min. 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.7713893  -0.0204082  0.0059170          10.2935       < 2.2e-16
    ## 2 (50)  0.4233853  -0.0204082  0.0035786           7.4187       9.464e-13
    ## 3 (50)  0.1689998  -0.0204082  0.0033251           3.2847        0.005105
    ## 4 (50)  0.0946069  -0.0204082  0.0043537           1.7431        0.243947
    ## 5 (50) -0.0270228  -0.0204082  0.0072209          -0.0778        0.937954
    ## 6 (40) -0.1788946  -0.0256410  0.0087697          -1.6365        0.243947
    ## 7 (30) -0.5970969  -0.0344828  0.0111079          -5.3382       5.633e-07
    ## 8 (20) -0.7668416  -0.0526316  0.0148195          -5.8669       3.108e-08
    ## 9 (10) -0.4866662  -0.1111111  0.0172835          -2.8567        0.017125
    ##           
    ## 1 (50) ***
    ## 2 (50) ***
    ## 3 (50) ** 
    ## 4 (50)    
    ## 5 (50)    
    ## 6 (40)    
    ## 7 (30) ***
    ## 8 (20) ***
    ## 9 (10) *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $pH
    ## Spatial correlogram for pH 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.7206808  -0.0204082  0.0060009           9.5667       < 2.2e-16
    ## 2 (50)  0.5123785  -0.0204082  0.0036285           8.8448       < 2.2e-16
    ## 3 (50)  0.2099793  -0.0204082  0.0033715           3.9678       0.0002902
    ## 4 (50) -0.3586281  -0.0204082  0.0044141          -5.0907       1.784e-06
    ## 5 (50) -0.6585951  -0.0204082  0.0073237          -7.4573       6.182e-13
    ## 6 (40) -0.6618167  -0.0256410  0.0089300          -6.7321       1.003e-10
    ## 7 (30) -0.3727208  -0.0344828  0.0113908          -3.1692       0.0045864
    ## 8 (20) -0.2124804  -0.0526316  0.0154419          -1.2864       0.3966415
    ## 9 (10) -0.1056156  -0.1111111  0.0193348           0.0395       0.9684743
    ##           
    ## 1 (50) ***
    ## 2 (50) ***
    ## 3 (50) ***
    ## 4 (50) ***
    ## 5 (50) ***
    ## 6 (40) ***
    ## 7 (30) ** 
    ## 8 (20)    
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $elevacion_media
    ## Spatial correlogram for elevacion_media 
    ## method: Moran's I
    ##           estimate expectation    variance standard deviate
    ## 1 (50)  0.74112321 -0.02040816  0.00601178           9.8217
    ## 2 (50)  0.20185740 -0.02040816  0.00363501           3.6865
    ## 3 (50) -0.26563804 -0.02040816  0.00337750          -4.2196
    ## 4 (50) -0.31713608 -0.02040816  0.00442198          -4.4622
    ## 5 (50) -0.01092311 -0.02040816  0.00733708           0.1107
    ## 6 (40) -0.01736264 -0.02564103  0.00895078           0.0875
    ## 7 (30) -0.00495651 -0.03448276  0.01142753           0.2762
    ## 8 (20) -0.00054486 -0.05263158  0.01552258           0.4181
    ## 9 (10)  0.00519599 -0.11111111  0.01960077           0.8307
    ##        Pr(I) two sided    
    ## 1 (50)       < 2.2e-16 ***
    ## 2 (50)       0.0013639 ** 
    ## 3 (50)       0.0001713 ***
    ## 4 (50)        6.49e-05 ***
    ## 5 (50)       1.0000000    
    ## 6 (40)       1.0000000    
    ## 7 (30)       1.0000000    
    ## 8 (20)       1.0000000    
    ## 9 (10)       1.0000000    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $pendiente_media
    ## Spatial correlogram for pendiente_media 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.2631057  -0.0204082  0.0060062           3.6583        0.002285
    ## 2 (50) -0.0135020  -0.0204082  0.0036317           0.1146        1.000000
    ## 3 (50) -0.0430490  -0.0204082  0.0033744          -0.3898        1.000000
    ## 4 (50) -0.1216382  -0.0204082  0.0044179          -1.5230        1.000000
    ## 5 (50) -0.0479517  -0.0204082  0.0073302          -0.3217        1.000000
    ## 6 (40) -0.0592347  -0.0256410  0.0089400          -0.3553        1.000000
    ## 7 (30) -0.0035707  -0.0344828  0.0114085           0.2894        1.000000
    ## 8 (20)  0.0432293  -0.0526316  0.0154808           0.7705        1.000000
    ## 9 (10)  0.0155156  -0.1111111  0.0194631           0.9077        1.000000
    ##          
    ## 1 (50) **
    ## 2 (50)   
    ## 3 (50)   
    ## 4 (50)   
    ## 5 (50)   
    ## 6 (40)   
    ## 7 (30)   
    ## 8 (20)   
    ## 9 (10)   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $orientacion_media
    ## Spatial correlogram for orientacion_media 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.4199983  -0.0204082  0.0057544           5.8057        5.77e-08
    ## 2 (50)  0.0561735  -0.0204082  0.0034818           1.2978          1.0000
    ## 3 (50)  0.0277699  -0.0204082  0.0032352           0.8470          1.0000
    ## 4 (50) -0.1550453  -0.0204082  0.0042367          -2.0685          0.3088
    ## 5 (50) -0.1218605  -0.0204082  0.0070215          -1.2107          1.0000
    ## 6 (40) -0.1064297  -0.0256410  0.0084592          -0.8784          1.0000
    ## 7 (30) -0.0436254  -0.0344828  0.0105597          -0.0890          1.0000
    ## 8 (20) -0.0547274  -0.0526316  0.0136137          -0.0180          1.0000
    ## 9 (10) -0.1021968  -0.1111111  0.0133092           0.0773          1.0000
    ##           
    ## 1 (50) ***
    ## 2 (50)    
    ## 3 (50)    
    ## 4 (50)    
    ## 5 (50)    
    ## 6 (40)    
    ## 7 (30)    
    ## 8 (20)    
    ## 9 (10)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $curvatura_perfil_media
    ## Spatial correlogram for curvatura_perfil_media 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50) -0.0935049  -0.0204082  0.0058877          -0.9526               1
    ## 2 (50) -0.0751714  -0.0204082  0.0035611          -0.9177               1
    ## 3 (50)  0.0217716  -0.0204082  0.0033089           0.7333               1
    ## 4 (50) -0.0201605  -0.0204082  0.0043326           0.0038               1
    ## 5 (50) -0.0181084  -0.0204082  0.0071849           0.0271               1
    ## 6 (40)  0.0045929  -0.0256410  0.0087137           0.3239               1
    ## 7 (30)  0.0046094  -0.0344828  0.0110090           0.3726               1
    ## 8 (20)  0.0522792  -0.0526316  0.0146018           0.8682               1
    ## 9 (10) -0.0337028  -0.1111111  0.0165661           0.6014               1
    ## 
    ## $curvatura_tangencial_media
    ## Spatial correlogram for curvatura_tangencial_media 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.0022552  -0.0204082  0.0057840           0.2980               1
    ## 2 (50) -0.0148778  -0.0204082  0.0034993           0.0935               1
    ## 3 (50) -0.0659636  -0.0204082  0.0032516          -0.7989               1
    ## 4 (50) -0.0079637  -0.0204082  0.0042580           0.1907               1
    ## 5 (50)  0.0074371  -0.0204082  0.0070577           0.3315               1
    ## 6 (40) -0.0997747  -0.0256410  0.0085156          -0.8034               1
    ## 7 (30) -0.0796818  -0.0344828  0.0106592          -0.4378               1
    ## 8 (20)  0.0281163  -0.0526316  0.0138325           0.6866               1
    ## 9 (10) -0.0020874  -0.1111111  0.0140305           0.9204               1
    ## 
    ## $abundancia_global
    ## Spatial correlogram for abundancia_global 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.1534056  -0.0204082  0.0058955           2.2637          0.2123
    ## 2 (50) -0.0129711  -0.0204082  0.0035658           0.1245          1.0000
    ## 3 (50) -0.0065438  -0.0204082  0.0033132           0.2409          1.0000
    ## 4 (50) -0.1558889  -0.0204082  0.0043382          -2.0569          0.3175
    ## 5 (50) -0.0483422  -0.0204082  0.0071945          -0.3293          1.0000
    ## 6 (40) -0.0185837  -0.0256410  0.0087286           0.0755          1.0000
    ## 7 (30) -0.0197171  -0.0344828  0.0110353           0.1406          1.0000
    ## 8 (20)  0.0271284  -0.0526316  0.0146598           0.6588          1.0000
    ## 9 (10)  0.0207136  -0.1111111  0.0167570           1.0184          1.0000
    ## 
    ## $riqueza_global
    ## Spatial correlogram for riqueza_global 
    ## method: Moran's I
    ##          estimate expectation   variance standard deviate Pr(I) two sided
    ## 1 (50)  0.1746908  -0.0204082  0.0057474           2.5735         0.09062
    ## 2 (50) -0.0158156  -0.0204082  0.0034776           0.0779         1.00000
    ## 3 (50) -0.1550748  -0.0204082  0.0032314          -2.3690         0.14269
    ## 4 (50) -0.1117792  -0.0204082  0.0042317          -1.4046         1.00000
    ## 5 (50)  0.0075414  -0.0204082  0.0070129           0.3338         1.00000
    ## 6 (40) -0.0531788  -0.0256410  0.0084458          -0.2996         1.00000
    ## 7 (30) -0.0389949  -0.0344828  0.0105361          -0.0440         1.00000
    ## 8 (20)  0.0081932  -0.0526316  0.0135616           0.5223         1.00000
    ## 9 (10)  0.0378480  -0.1111111  0.0131376           1.2996         1.00000
    ##         
    ## 1 (50) .
    ## 2 (50)  
    ## 3 (50)  
    ## 4 (50)  
    ## 5 (50)  
    ## 6 (40)  
    ## 7 (30)  
    ## 8 (20)  
    ## 9 (10)  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
dim_panel <- rev(n2mfrow(ncol(bci_env_grid_num)))
```

``` r
par(mfrow = dim_panel)
suppressWarnings(invisible(lapply(auto_amb, function(x) plot(x, main = x$var))))
```

![](ee_ecologia_espacial_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

### Correlograma Mantel de matriz datos comunidad sin tendencia (residuos)

``` r
mi_fam_sin_tendencia <- resid(
  lm(as.matrix(mi_fam_hel) ~ .,
     data = bci_xy))
mi_fam_sin_tendencia_d <- dist(mi_fam_sin_tendencia)
(mi_fam_correlograma <- mantel.correlog(
  mi_fam_sin_tendencia_d,
  XY = bci_xy,
  nperm = 999))
```

    ## 
    ## Mantel Correlogram Analysis
    ## 
    ## Call:
    ##  
    ## mantel.correlog(D.eco = mi_fam_sin_tendencia_d, XY = bci_xy,      nperm = 999) 
    ## 
    ##         class.index     n.dist Mantel.cor Pr(Mantel) Pr(corrected)   
    ## D.cl.1   136.870241 144.000000   0.064513      0.002         0.002 **
    ## D.cl.2   210.610723 376.000000   0.038638      0.091         0.091 . 
    ## D.cl.3   284.351204 390.000000  -0.052374      0.026         0.052 . 
    ## D.cl.4   358.091686 148.000000  -0.042984      0.023         0.069 . 
    ## D.cl.5   431.832168 372.000000  -0.059022      0.023         0.092 . 
    ## D.cl.6   505.572649 266.000000  -0.052373      0.004         0.020 * 
    ## D.cl.7   579.313131 168.000000         NA         NA            NA   
    ## D.cl.8   653.053613 100.000000         NA         NA            NA   
    ## D.cl.9   726.794094 154.000000         NA         NA            NA   
    ## D.cl.10  800.534576  88.000000         NA         NA            NA   
    ## D.cl.11  874.275058  50.000000         NA         NA            NA   
    ## D.cl.12  948.015539  24.000000         NA         NA            NA   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(mi_fam_correlograma)
```

![](ee_ecologia_espacial_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

### Determinación de autocorrelación global de residuos por medio de prueba de permutación del I de Moran

``` r
autocor_global_residuos <- sapply(
  dimnames(mi_fam_sin_tendencia)[[2]],
  function(x)
    moran.mc(
      x = mi_fam_sin_tendencia[,x],
      listw = pesos_b,
      zero.policy = T,
      nsim = 9999),
    simplify = F)

bci_env_grid_num_sf <- bci_env_grid %>%
  select_if(is.numeric) %>% 
  select(-id, -UTM.EW, -UTM.NS)
```

### Determinación de autocorrelación local de variables ambientales por medio de prueba de I de Moran

``` r
lisamaps_amb <- sapply(grep('geometry', names(bci_env_grid_num_sf), invert = T, value = T),
                   function(x) {
                     m <- lisamap(objesp = bci_env_grid_num_sf[x],
                                  var = x,
                                  pesos = pesos_b,
                                  tituloleyenda = 'Significancia ("x-y", léase como "x" rodeado de "y")',
                                  leyenda = F,
                                  anchuratitulo = 50,
                                  tamanotitulo = 10,
                                  fuentedatos = '\nhttp://ctfs.si.edu/webatlas/datasets/bci/',
                                  titulomapa = paste0('Clusters LISA de "', x, '"'))
                     return(m$grafico)
                   }, simplify = F
)
lisamaps_amb$leyenda <- gtable_filter(ggplot_gtable(ggplot_build(lisamaps_amb[[1]] + theme(legend.position="bottom"))), "guide-box")
grid.arrange(do.call('arrangeGrob', c(lisamaps_amb[1:12], nrow = 3)), lisamaps_amb$leyenda, heights=c(1.1, 0.1), nrow = 2)
```

![](ee_ecologia_espacial_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
grid.arrange(do.call('arrangeGrob', c(lisamaps_amb[13:22], nrow = 3)), lisamaps_amb$leyenda, heights=c(1.1, 0.1), nrow = 2)
```

![](ee_ecologia_espacial_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

``` r
grid.arrange(do.call('arrangeGrob', c(lisamaps_amb[23:31], nrow = 3)), lisamaps_amb$leyenda, heights=c(1.1, 0.1), nrow = 2)
```

![](ee_ecologia_espacial_files/figure-gfm/unnamed-chunk-16-3.png)<!-- -->

### Determinación de autocorrelación local de abundancias transformadas (Hellinger) por medio de prueba de I de Moran

``` r
mi_fam_hel_sf <- bci_env_grid %>% select %>% bind_cols(mi_fam_hel)
lisamaps_mifam <- sapply(
  grep('geometry', names(mi_fam_hel_sf), invert = T, value = T),
                   function(x) {
                     m <- lisamap(objesp = mi_fam_hel_sf[x],
                                  var = x,
                                  pesos = pesos_b,
                                  tituloleyenda = 'Significancia ("x-y", léase como "x" rodeado de "y")',
                                  leyenda = F,
                                  anchuratitulo = 50,
                                  tamanotitulo = 10,
                                  fuentedatos = '\nhttp://ctfs.si.edu/webatlas/datasets/bci/',
                                  titulomapa = paste0('Clusters LISA de "', x, '"'))
                     # dev.new();print(m$grafico)
                     return(m$grafico)
                   }, simplify = F
)
lisamaps_mifam$leyenda <- gtable_filter(ggplot_gtable(ggplot_build(lisamaps_mifam[[1]] + theme(legend.position="bottom"))), "guide-box")
grid.arrange(do.call('arrangeGrob', c(lisamaps_mifam[1:8], nrow = 3)), lisamaps_mifam$leyenda, heights=c(1.1, 0.1), nrow = 2)
```

![](ee_ecologia_espacial_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
grid.arrange(do.call('arrangeGrob', c(lisamaps_mifam[9:16], nrow = 3)), lisamaps_mifam$leyenda, heights=c(1.1, 0.1), nrow = 2)
```

![](ee_ecologia_espacial_files/figure-gfm/unnamed-chunk-17-2.png)<!-- -->
