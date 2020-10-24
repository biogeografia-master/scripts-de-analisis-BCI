Análisis exploratorio de datos. Colección tidyverse
================
JR
18 de octubre, 2020

# ¿Qué es tidyverse?

Es una colección de paquetes con los que podrás importar, transformar,
visualizar, modelar y presentar datos. La colección se compone de 8
paquetes, de los cuales verás sobre todo 4: `ggplot2`, `dplyr`, `tidyr`
y `readr`.

Todos estos paquetes comparten estructuras comunes. Una de las
herramientas que incorpora la colección es la pipa `%>%` (**SHORTCUT:
`CTRL+SHIFT+M`**), la cual importa desde el paquete `magrittr`. Usarás
la pipa para construir “tuberías” de procesamiento sin necesidad de
crear objetos intermedios. En una tubería, puedes interpretar la pipa
como **“luego”**, y verás más adelante por qué. La función principal de
la pipa (tiene muchas, pero esta es la más importante) es pasar el
resultado del objeto a su izquierda como primer argumento de la función
a su derecha. El siguiente ejemplo explica su uso:

`objeto1 %>% funcion1()` es equivalente a `funcion1(argumento1 =
objeto1)`

> La idea del *pipe* pertenece a la tradición de sistemas Unix y, en
> origen, su función era comunicar distintos procesos, usando la salida
> estándar de uno (*stdout*) como entrada estándar (*stdin*) del
> siguiente.

Su ventaja radica en que, si necesitaras continuar procesando los datos,
no tendrás que anidar ni crear objetos intermedios. En el siguiente
ejemplo, asigno el resultado de una cadena al objeto nombrado
`resultado`:

`resultado <- objeto1 %>% funcion1() %>% funcion2() %>% funcion3()`

Puedes leer lo anterior como *"objeto1 pasa como primer argumento de
funcion1, **luego** el resultado de funcion1 pasa como primer argumento
de funcion2, **luego** el resultado de funcion2 pasa como primer
argumento de funcion3*.

Para replicar esta operación sin la pipa, se necesitaría algo tal que
ésto:

`resultado <- funcion3(funcion2(funcion1(objeto1)))`

O alternativamente, crear objetos intermedios:

`tmp1 <- funcion1(objeto1)` `tmp2 <- funcion2(tmp1)` `resultado <-
funcion3(tmp2)`

Notarás que la tubería es más limpia que estas dos últimas opciones. La
tubería puedes leerla de forma encandenada, a diferencia del estilo
anidado y de creación de objetos intermedios, que añade una cierta
complejidad de lectura para el usuario/a, sobre todo para personas sin
conocimientos de programación. Precisamente por esta razón fue que
decidí introducir la colección tidyverse, para así mostrarte algunas
ideas que podrás aplicar a tus datos y, en principio, para facilitarte
la vida (“*no me ayude’ tali*”). Ahora bien, si decides programar en R
más adelante, deberás aprender las capacidades de programación
orientada a objetos y programación funcional de
    R.

¡Comencemos\!

## Paquetes

``` r
library(tidyverse)
```

    ## ── Attaching packages ──────────────────────────────────── tidyverse 1.2.1 ──

    ## ✓ ggplot2 3.3.2     ✓ purrr   0.3.4
    ## ✓ tibble  3.0.3     ✓ dplyr   0.8.3
    ## ✓ tidyr   1.0.0     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.4.0

    ## ── Conflicts ─────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(sf)
```

    ## Linking to GEOS 3.6.2, GDAL 2.2.3, PROJ 4.9.3

> `sf` te ayudará a leer el objeto `bci_env_grid` como un *simple
> feature*, el cual se encuentra dentro del archivo
> `biodata/matriz_ambiental.Rdata`. Esto extenderá las capacidades
> espaciales del objeto.

## Cargar datos

``` r
load('biodata/matriz_ambiental.Rdata')
load('biodata/Apocynaceae-Meliaceae-Sapotaceae.Rdata')
```

## Paquete `dplyr`

Te servirá para manipular datos mediante verbos. Los verbos de `dplyr`
que conocerás son (hay muchos otros): `select()`, `filter()`,
`arrange()`, `mutate()`, `group_by()` y `summarise()`.

### Verbo `select`

Comúnmente, necesitas seleccionar una o varias columnas de una tabla.
Para esto existe el verbo `select`. Te muestro un ejemplo aplicado a la
matriz de comunidad (que por ahora la verás como `simple feature`),
seleccionando las columnas `id` (número identificador de quadrat) y `pH`
(pH del suelo):

``` r
bci_env_grid %>%
  select(id, pH)
```

    ## Simple feature collection with 50 features and 2 fields
    ## geometry type:  POLYGON
    ## dimension:      XY
    ## bbox:           xmin: 625704 ymin: 1011519 xmax: 626704 ymax: 1012019
    ## epsg (SRID):    32617
    ## proj4string:    +proj=utm +zone=17 +datum=WGS84 +units=m +no_defs
    ## First 10 features:
    ##    id      pH                       geometry
    ## 1   1 4.32432 POLYGON ((625704 1011519, 6...
    ## 2   2 4.37548 POLYGON ((625704 1011619, 6...
    ## 3   3 4.34700 POLYGON ((625704 1011719, 6...
    ## 4   4 4.46112 POLYGON ((625704 1011819, 6...
    ## 5   5 4.40128 POLYGON ((625704 1011919, 6...
    ## 6   6 4.57252 POLYGON ((625804 1011519, 6...
    ## 7   7 4.55972 POLYGON ((625804 1011619, 6...
    ## 8   8 4.41168 POLYGON ((625804 1011719, 6...
    ## 9   9 4.53336 POLYGON ((625804 1011819, 6...
    ## 10 10 4.55500 POLYGON ((625804 1011919, 6...

> Importante: el objeto `bci_env_grid` permanece intacto, a menos que se
> use dicho nombre para reasignarlo a otro objeto. Mientras no se use el
> asignador `<-`, sólo verás que manipulo y visualizo copias del objeto
> original. Fíjate en la clase del objeto `bci_env_grid`:

``` r
bci_env_grid %>%
  class
```

    ## [1] "sf"         "data.frame"

El objeto `bci_env_grid` es a la vez de clase `sf` (*simple feature*) y
`data.frame`, es decir, es tanto tabla como objeto espacial, por lo que
se puede representar en un mapa. Este objeto no pierde la clase `sf`,
por lo que verás que aparece información geométrica y geoespacial en el
encabezado, y luego un extracto de la tabla de datos (como máximo, las
10 primeras filas). Para convertirlo a un simple `data.frame`, hay que
“tumbar” su geometría con `st_drop_geometry`:

``` r
bci_env_grid %>%
  select(id, pH) %>%
  st_drop_geometry
```

    ##    id      pH
    ## 1   1 4.32432
    ## 2   2 4.37548
    ## 3   3 4.34700
    ## 4   4 4.46112
    ## 5   5 4.40128
    ## 6   6 4.57252
    ## 7   7 4.55972
    ## 8   8 4.41168
    ## 9   9 4.53336
    ## 10 10 4.55500
    ## 11 11 4.71792
    ## 12 12 4.29640
    ## 13 13 4.12084
    ## 14 14 4.12820
    ## 15 15 4.15760
    ## 16 16 4.50452
    ## 17 17 4.46092
    ## 18 18 4.05204
    ## 19 19 4.24652
    ## 20 20 4.33736
    ## 21 21 4.77492
    ## 22 22 4.67260
    ## 23 23 4.29220
    ## 24 24 4.69712
    ## 25 25 4.74796
    ## 26 26 4.96144
    ## 27 27 4.97480
    ## 28 28 4.84068
    ## 29 29 4.75024
    ## 30 30 4.75976
    ## 31 31 5.00428
    ## 32 32 5.15908
    ## 33 33 5.05812
    ## 34 34 5.05132
    ## 35 35 4.85808
    ## 36 36 4.93364
    ## 37 37 4.93540
    ## 38 38 4.88152
    ## 39 39 4.93488
    ## 40 40 4.77104
    ## 41 41 4.89252
    ## 42 42 4.82900
    ## 43 43 4.82804
    ## 44 44 4.94864
    ## 45 45 5.02364
    ## 46 46 4.68412
    ## 47 47 4.52992
    ## 48 48 4.87296
    ## 49 49 5.00388
    ## 50 50 5.02052

Fíjate ahora en la clase de `bci_env_grid %>% select(id, pH) %>%
st_drop_geometry`, que en este caso es sólo `data.frame`:

``` r
bci_env_grid %>%
  select(id, pH) %>%
  st_drop_geometry %>%
  class
```

    ## [1] "data.frame"

> Al introducir un `<enter>` después de la pipa, el código puede
> continuar en la línea siguiente. Esto se hace para evitar que la línea
> de código sea legible sin necesidad de desplazarse hacia la derecha
> (es aconsejable no superar los 80 caracteres en una misma línea de
> código, según recomiendan en la [Google’s R Style
> Guide](https://google.github.io/styleguide/Rguide.html)). Como
> convención, escribiré un `<enter>` después de cada operador pipa.

Seleccionaré, y a la vez renombraré, dos columnas con `select`
(recuerda: no estoy modificando el objeto original, simplemente trabajo
en copias no asignadas). De paso, sólo mostraré las 6 primeras filas al
aplicar `head` al final de la tubería (no sólo se admiten verbos
`dplyr`, cualquier función de R puede entrar en la tubería):

``` r
bci_env_grid %>%
  select(id_de_quadrat = id, pH_del_suelo = pH) %>%
  st_drop_geometry %>%
  head
```

    ##   id_de_quadrat pH_del_suelo
    ## 1             1      4.32432
    ## 2             2      4.37548
    ## 3             3      4.34700
    ## 4             4      4.46112
    ## 5             5      4.40128
    ## 6             6      4.57252

Ahora mostraré sólo los elementos con `pH` mayor que 5, usando el verbo
`filter`

``` r
bci_env_grid %>%
  select(id, pH) %>%
  st_drop_geometry %>%
  filter(pH>5)
```

    ##   id      pH
    ## 1 31 5.00428
    ## 2 32 5.15908
    ## 3 33 5.05812
    ## 4 34 5.05132
    ## 5 45 5.02364
    ## 6 49 5.00388
    ## 7 50 5.02052

O filtro por aquellos con `id` 31 y 50:

``` r
bci_env_grid %>%
  select(id, pH) %>%
  st_drop_geometry %>%
  filter(id == c(31, 50))
```

    ##   id      pH
    ## 1 31 5.00428
    ## 2 50 5.02052

Pruebo también con la matriz de comunidad. Por ejemplo, introduzco en la
tubería la función `colSums`, que devuelve un vector cuyos elementos
están nombrados (tienen un atributo, en este caso, el nombre de
especie), donde cada elemento representa la abundancia por especie.

``` r
mc_apcyn_melic_saptc %>%
  colSums
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

Y también obtengo la abundancia por quadrat.

``` r
mc_apcyn_melic_saptc %>%
  rowSums
```

    ##   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18 
    ## 309 313 215 237 188 193 253 215 180 214 266 231 306 278 315 349 217 385 
    ##  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36 
    ## 201 429 282 234 491 219 551 371 414 311 254 459 432 426 610 581 531 350 
    ##  37  38  39  40  41  42  43  44  45  46  47  48  49  50 
    ## 340 359 555 593 361 298 430 684 609 167 257 724 745 494

Uso a continuación el verbo `arrange` para mostrar los registros de la
matriz ambiental ordenados ascendentemente por pH.

``` r
bci_env_grid %>%
  select(id, pH) %>%
  st_drop_geometry %>%
  arrange(pH)
```

    ##    id      pH
    ## 1  18 4.05204
    ## 2  13 4.12084
    ## 3  14 4.12820
    ## 4  15 4.15760
    ## 5  19 4.24652
    ## 6  23 4.29220
    ## 7  12 4.29640
    ## 8   1 4.32432
    ## 9  20 4.33736
    ## 10  3 4.34700
    ## 11  2 4.37548
    ## 12  5 4.40128
    ## 13  8 4.41168
    ## 14 17 4.46092
    ## 15  4 4.46112
    ## 16 16 4.50452
    ## 17 47 4.52992
    ## 18  9 4.53336
    ## 19 10 4.55500
    ## 20  7 4.55972
    ## 21  6 4.57252
    ## 22 22 4.67260
    ## 23 46 4.68412
    ## 24 24 4.69712
    ## 25 11 4.71792
    ## 26 25 4.74796
    ## 27 29 4.75024
    ## 28 30 4.75976
    ## 29 40 4.77104
    ## 30 21 4.77492
    ## 31 43 4.82804
    ## 32 42 4.82900
    ## 33 28 4.84068
    ## 34 35 4.85808
    ## 35 48 4.87296
    ## 36 38 4.88152
    ## 37 41 4.89252
    ## 38 36 4.93364
    ## 39 39 4.93488
    ## 40 37 4.93540
    ## 41 44 4.94864
    ## 42 26 4.96144
    ## 43 27 4.97480
    ## 44 49 5.00388
    ## 45 31 5.00428
    ## 46 50 5.02052
    ## 47 45 5.02364
    ## 48 34 5.05132
    ## 49 33 5.05812
    ## 50 32 5.15908

Ahora usaré `arrange` para mostrar los registros de la matriz ambiental
ordenados DESCendentemente por pH.

``` r
bci_env_grid %>%
  select(id, pH) %>%
  st_drop_geometry %>%
  arrange(desc(pH))
```

    ##    id      pH
    ## 1  32 5.15908
    ## 2  33 5.05812
    ## 3  34 5.05132
    ## 4  45 5.02364
    ## 5  50 5.02052
    ## 6  31 5.00428
    ## 7  49 5.00388
    ## 8  27 4.97480
    ## 9  26 4.96144
    ## 10 44 4.94864
    ## 11 37 4.93540
    ## 12 39 4.93488
    ## 13 36 4.93364
    ## 14 41 4.89252
    ## 15 38 4.88152
    ## 16 48 4.87296
    ## 17 35 4.85808
    ## 18 28 4.84068
    ## 19 42 4.82900
    ## 20 43 4.82804
    ## 21 21 4.77492
    ## 22 40 4.77104
    ## 23 30 4.75976
    ## 24 29 4.75024
    ## 25 25 4.74796
    ## 26 11 4.71792
    ## 27 24 4.69712
    ## 28 46 4.68412
    ## 29 22 4.67260
    ## 30  6 4.57252
    ## 31  7 4.55972
    ## 32 10 4.55500
    ## 33  9 4.53336
    ## 34 47 4.52992
    ## 35 16 4.50452
    ## 36  4 4.46112
    ## 37 17 4.46092
    ## 38  8 4.41168
    ## 39  5 4.40128
    ## 40  2 4.37548
    ## 41  3 4.34700
    ## 42 20 4.33736
    ## 43  1 4.32432
    ## 44 12 4.29640
    ## 45 23 4.29220
    ## 46 19 4.24652
    ## 47 15 4.15760
    ## 48 14 4.12820
    ## 49 13 4.12084
    ## 50 18 4.05204

Usaré el verbo `mutate` para crear una nueva columna. Por ejemplo, creo
una columna que contenga `habitat` y `quebrada` separadas por una coma:

``` r
bci_env_grid %>%
  st_drop_geometry %>%
  select(habitat, quebrada) %>% 
  mutate(habitat_quebrada = paste(habitat, quebrada, sep = ', '))
```

    ##     habitat quebrada habitat_quebrada
    ## 1  OldSlope      Yes    OldSlope, Yes
    ## 2    OldLow      Yes      OldLow, Yes
    ## 3    OldLow       No       OldLow, No
    ## 4    OldLow       No       OldLow, No
    ## 5  OldSlope       No     OldSlope, No
    ## 6    OldLow       No       OldLow, No
    ## 7    OldLow      Yes      OldLow, Yes
    ## 8    OldLow      Yes      OldLow, Yes
    ## 9    OldLow       No       OldLow, No
    ## 10   OldLow       No       OldLow, No
    ## 11   OldLow       No       OldLow, No
    ## 12   OldLow       No       OldLow, No
    ## 13   OldLow      Yes      OldLow, Yes
    ## 14   OldLow       No       OldLow, No
    ## 15   OldLow       No       OldLow, No
    ## 16 OldSlope       No     OldSlope, No
    ## 17   OldLow       No       OldLow, No
    ## 18    Swamp       No        Swamp, No
    ## 19   OldLow       No       OldLow, No
    ## 20   OldLow       No       OldLow, No
    ## 21 OldSlope       No     OldSlope, No
    ## 22   OldLow       No       OldLow, No
    ## 23    Swamp       No        Swamp, No
    ## 24   OldLow       No       OldLow, No
    ## 25   OldLow       No       OldLow, No
    ## 26 OldSlope       No     OldSlope, No
    ## 27   OldLow       No       OldLow, No
    ## 28   OldLow       No       OldLow, No
    ## 29  OldHigh       No      OldHigh, No
    ## 30    Young       No        Young, No
    ## 31   OldLow       No       OldLow, No
    ## 32  OldHigh       No      OldHigh, No
    ## 33  OldHigh       No      OldHigh, No
    ## 34  OldHigh       No      OldHigh, No
    ## 35    Young       No        Young, No
    ## 36 OldSlope       No     OldSlope, No
    ## 37  OldHigh       No      OldHigh, No
    ## 38  OldHigh       No      OldHigh, No
    ## 39  OldHigh       No      OldHigh, No
    ## 40  OldHigh       No      OldHigh, No
    ## 41 OldSlope       No     OldSlope, No
    ## 42 OldSlope       No     OldSlope, No
    ## 43 OldSlope       No     OldSlope, No
    ## 44 OldSlope       No     OldSlope, No
    ## 45 OldSlope      Yes    OldSlope, Yes
    ## 46   OldLow       No       OldLow, No
    ## 47   OldLow       No       OldLow, No
    ## 48   OldLow       No       OldLow, No
    ## 49   OldLow       No       OldLow, No
    ## 50 OldSlope      Yes    OldSlope, Yes

Ahora `mutate`, pero con números: creo una columna de área de cada
cuadro (necesitas también la función `st_area`, del paquete `sf`):

``` r
bci_env_grid %>%
  mutate(area = st_area(geometry)) %>%
  select(id, area) %>%
  st_drop_geometry %>%
  head
```

    ##   id        area
    ## 1  1 10000 [m^2]
    ## 2  2 10000 [m^2]
    ## 3  3 10000 [m^2]
    ## 4  4 10000 [m^2]
    ## 5  5 10000 [m^2]
    ## 6  6 10000 [m^2]

…y ahora más complejo: obtengo la densidad de individuos por metro
cuadrado, ordenados descendentemente por dicha densidad, y conservando
sólo los 6 registros con mayores densidades.

``` r
bci_env_grid %>%
  mutate(area = st_area(geometry), densidad_indiv = abundancia_global/area) %>% 
  select(id, densidad_indiv) %>% 
  st_drop_geometry %>% 
  arrange(desc(densidad_indiv)) %>% 
  head
```

    ##   id densidad_indiv
    ## 1 20 0.5014 [1/m^2]
    ## 2  3 0.4611 [1/m^2]
    ## 3 48 0.4597 [1/m^2]
    ## 4 15 0.4559 [1/m^2]
    ## 5  5 0.4549 [1/m^2]
    ## 6 43 0.4531 [1/m^2]

Finalmente, te muestro los verbos `group_by` y `summarise`, los cuales
son útiles para producir resúmenes por grupos. Agruparé la matriz
ambiental por la columna `habitat`, dejando sólo las variables numericas
que hagan sentido (por ejemplo, excluyo `id`, `UTM.EW`, `UTM.NS`), y lo
asignaré a `agrupado_por_habitat` para luego reutilizarlo:

``` r
agrupado_por_habitat <- bci_env_grid %>%
  st_drop_geometry %>%
  group_by(habitat) %>% 
  select_if(is.numeric) %>%
  select(-id, -UTM.EW, -UTM.NS)
agrupado_por_habitat
```

    ## # A tibble: 50 x 32
    ## # Groups:   habitat [5]
    ##    habitat heterogeneidad_… geomorf_llanura… geomorf_pico_pct
    ##  * <fct>              <dbl>            <dbl>            <dbl>
    ##  1 OldSlo…           0.627             10.0              0   
    ##  2 OldLow            0.394             34.8              0   
    ##  3 OldLow            0                  0                0   
    ##  4 OldLow            0                  0                0   
    ##  5 OldSlo…           0.461              2.58             0   
    ##  6 OldLow            0.0768             0                0.17
    ##  7 OldLow            0.381              0                0.53
    ##  8 OldLow            0.211              0                0   
    ##  9 OldLow            0                  0                0   
    ## 10 OldLow            0                  1.03             0   
    ## # … with 40 more rows, and 28 more variables:
    ## #   geomorf_interfluvio_pct <dbl>, geomorf_hombrera_pct <dbl>,
    ## #   `geomorf_espolón/gajo_pct` <dbl>, geomorf_vertiente_pct <dbl>,
    ## #   geomorf_vaguada_pct <dbl>, geomorf_piedemonte_pct <dbl>,
    ## #   geomorf_valle_pct <dbl>, geomorf_sima_pct <dbl>, Al <dbl>, B <dbl>,
    ## #   Ca <dbl>, Cu <dbl>, Fe <dbl>, K <dbl>, Mg <dbl>, Mn <dbl>, P <dbl>,
    ## #   Zn <dbl>, N <dbl>, N.min. <dbl>, pH <dbl>, elevacion_media <dbl>,
    ## #   pendiente_media <dbl>, orientacion_media <dbl>,
    ## #   curvatura_perfil_media <dbl>, curvatura_tangencial_media <dbl>,
    ## #   abundancia_global <dbl>, riqueza_global <int>

Observa el encabezado: el objeto es `A tibble: 50 x 32` y hay 5 grupos
(`Groups: habitat [5]`). Calculo cuántos elementos (filas) hay por
grupo:

``` r
agrupado_por_habitat %>% summarise(n = n())
```

    ## # A tibble: 5 x 2
    ##   habitat      n
    ##   <fct>    <int>
    ## 1 OldHigh      8
    ## 2 OldLow      26
    ## 3 OldSlope    12
    ## 4 Swamp        2
    ## 5 Young        2

…y también algunos estadísticos de las columnas `pH`,
`abundancia_global` y `riqueza_global` por ejemplo:

``` r
agrupado_por_habitat %>%
  summarise(
    n = n(),
    media_pH = mean(pH),
    media_abundancia = mean(abundancia_global),
    media_riqueza = mean(riqueza_global)
  )
```

    ## # A tibble: 5 x 5
    ##   habitat      n media_pH media_abundancia media_riqueza
    ##   <fct>    <int>    <dbl>            <dbl>         <dbl>
    ## 1 OldHigh      8     4.94            3969.          157.
    ## 2 OldLow      26     4.55            4182.          168.
    ## 3 OldSlope    12     4.79            4036.          167.
    ## 4 Swamp        2     4.17            3832.          188.
    ## 5 Young        2     4.81            3882           160

…o la media de todas las variables numéricas

``` r
agrupado_por_habitat %>%
  summarise_all(mean)
```

    ## # A tibble: 5 x 32
    ##   habitat heterogeneidad_… geomorf_llanura… geomorf_pico_pct
    ##   <fct>              <dbl>            <dbl>            <dbl>
    ## 1 OldHigh            0.305            15.5            0     
    ## 2 OldLow             0.235            14.1            0.0308
    ## 3 OldSlo…            0.396             6.90           0     
    ## 4 Swamp              0.642             1.7            0     
    ## 5 Young              0.474            32.4            0     
    ## # … with 28 more variables: geomorf_interfluvio_pct <dbl>,
    ## #   geomorf_hombrera_pct <dbl>, `geomorf_espolón/gajo_pct` <dbl>,
    ## #   geomorf_vertiente_pct <dbl>, geomorf_vaguada_pct <dbl>,
    ## #   geomorf_piedemonte_pct <dbl>, geomorf_valle_pct <dbl>,
    ## #   geomorf_sima_pct <dbl>, Al <dbl>, B <dbl>, Ca <dbl>, Cu <dbl>,
    ## #   Fe <dbl>, K <dbl>, Mg <dbl>, Mn <dbl>, P <dbl>, Zn <dbl>, N <dbl>,
    ## #   N.min. <dbl>, pH <dbl>, elevacion_media <dbl>, pendiente_media <dbl>,
    ## #   orientacion_media <dbl>, curvatura_perfil_media <dbl>,
    ## #   curvatura_tangencial_media <dbl>, abundancia_global <dbl>,
    ## #   riqueza_global <dbl>

…no caben, mejor por partes

``` r
agrupado_por_habitat %>%
  summarise_all(mean) %>% 
  select(1:6) %>% 
  print(width=300)
```

    ## # A tibble: 5 x 6
    ##   habitat  heterogeneidad_ambiental geomorf_llanura_pct geomorf_pico_pct
    ##   <fct>                       <dbl>               <dbl>            <dbl>
    ## 1 OldHigh                     0.305               15.5            0     
    ## 2 OldLow                      0.235               14.1            0.0308
    ## 3 OldSlope                    0.396                6.90           0     
    ## 4 Swamp                       0.642                1.7            0     
    ## 5 Young                       0.474               32.4            0     
    ##   geomorf_interfluvio_pct geomorf_hombrera_pct
    ##                     <dbl>                <dbl>
    ## 1                   0.53                  3.21
    ## 2                   0.962                 2.85
    ## 3                   0.48                  2.21
    ## 4                   0                     2.4 
    ## 5                   0                     3.42

``` r
agrupado_por_habitat %>%
  summarise_all(mean) %>% 
  select(1,7:12) %>% 
  print(width=300)
```

    ## # A tibble: 5 x 7
    ##   habitat  `geomorf_espolón/gajo_pct` geomorf_vertiente_pct
    ##   <fct>                         <dbl>                 <dbl>
    ## 1 OldHigh                       2.80                   72.8
    ## 2 OldLow                        4.10                   70.9
    ## 3 OldSlope                      5.71                   80.9
    ## 4 Swamp                         1.32                   92.1
    ## 5 Young                         0.105                  56.8
    ##   geomorf_vaguada_pct geomorf_piedemonte_pct geomorf_valle_pct
    ##                 <dbl>                  <dbl>             <dbl>
    ## 1               1.75                   2.80              0.61 
    ## 2               3.26                   1.81              1.90 
    ## 3               2.85                   0.557             0.422
    ## 4               1.78                   0.355             0.38 
    ## 5               0.355                  6.59              0.4  
    ##   geomorf_sima_pct
    ##              <dbl>
    ## 1          0      
    ## 2          0.00423
    ## 3          0      
    ## 4          0      
    ## 5          0

``` r
agrupado_por_habitat %>%
  summarise_all(mean) %>% 
  select(1, 13:25) %>% 
  print(width=300)
```

    ## # A tibble: 5 x 14
    ##   habitat     Al     B    Ca    Cu    Fe     K    Mg    Mn     P    Zn
    ##   <fct>    <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ## 1 OldHigh  1043. 1.18  1717.  7.03  164. 182.   295.  401.  4.85  6.25
    ## 2 OldLow   1046. 0.737 1558.  6.62  167. 147.   276.  350.  2.45  4.29
    ## 3 OldSlope  883. 1.25  2135.  8.50  211. 216.   360.  438.  2.11  6.73
    ## 4 Swamp    1205. 0.364 1016.  5.04  198.  98.3  215.  117.  4.85  2.81
    ## 5 Young    1036. 1.08  1481.  6.54  164. 146.   264.  305.  3.11  4.64
    ##       N N.min.    pH
    ##   <dbl>  <dbl> <dbl>
    ## 1  34.2   18.2  4.94
    ## 2  22.2   15.5  4.55
    ## 3  27.8   23.6  4.79
    ## 4  20.2   16.0  4.17
    ## 5  34.7   14.2  4.81

``` r
agrupado_por_habitat %>%
  summarise_all(mean) %>% 
  select(1, 26:32) %>%
  print(width=300)
```

    ## # A tibble: 5 x 8
    ##   habitat  elevacion_media pendiente_media orientacion_media
    ##   <fct>              <dbl>           <dbl>             <dbl>
    ## 1 OldHigh             153.            6.46              237.
    ## 2 OldLow              142.            6.88              228.
    ## 3 OldSlope            140.            8.31              263.
    ## 4 Swamp               146.           12.1               278.
    ## 5 Young               161.            3.66              186.
    ##   curvatura_perfil_media curvatura_tangencial_media abundancia_global
    ##                    <dbl>                      <dbl>             <dbl>
    ## 1              0.000382                   0.000245              3969.
    ## 2              0.0000968                 -0.0000619             4182.
    ## 3              0.000848                   0.000277              4036.
    ## 4              0.000663                  -0.000187              3832.
    ## 5             -0.000192                  -0.0000480             3882 
    ##   riqueza_global
    ##            <dbl>
    ## 1           157.
    ## 2           168.
    ## 3           167.
    ## 4           188.
    ## 5           160

…y no sólo un estadístico, sino varios:

``` r
agrupado_por_habitat %>%
  summarise_all(
    list(
      media = mean,
      mediana = median,
      varianza = var,
      minimo = min,
      maximo = max
    )
  )
```

    ## # A tibble: 5 x 156
    ##   habitat heterogeneidad_… geomorf_llanura… geomorf_pico_pc…
    ##   <fct>              <dbl>            <dbl>            <dbl>
    ## 1 OldHigh            0.305            15.5            0     
    ## 2 OldLow             0.235            14.1            0.0308
    ## 3 OldSlo…            0.396             6.90           0     
    ## 4 Swamp              0.642             1.7            0     
    ## 5 Young              0.474            32.4            0     
    ## # … with 152 more variables: geomorf_interfluvio_pct_media <dbl>,
    ## #   geomorf_hombrera_pct_media <dbl>,
    ## #   `geomorf_espolón/gajo_pct_media` <dbl>,
    ## #   geomorf_vertiente_pct_media <dbl>, geomorf_vaguada_pct_media <dbl>,
    ## #   geomorf_piedemonte_pct_media <dbl>, geomorf_valle_pct_media <dbl>,
    ## #   geomorf_sima_pct_media <dbl>, Al_media <dbl>, B_media <dbl>,
    ## #   Ca_media <dbl>, Cu_media <dbl>, Fe_media <dbl>, K_media <dbl>,
    ## #   Mg_media <dbl>, Mn_media <dbl>, P_media <dbl>, Zn_media <dbl>,
    ## #   N_media <dbl>, N.min._media <dbl>, pH_media <dbl>,
    ## #   elevacion_media_media <dbl>, pendiente_media_media <dbl>,
    ## #   orientacion_media_media <dbl>, curvatura_perfil_media_media <dbl>,
    ## #   curvatura_tangencial_media_media <dbl>, abundancia_global_media <dbl>,
    ## #   riqueza_global_media <dbl>, heterogeneidad_ambiental_mediana <dbl>,
    ## #   geomorf_llanura_pct_mediana <dbl>, geomorf_pico_pct_mediana <dbl>,
    ## #   geomorf_interfluvio_pct_mediana <dbl>,
    ## #   geomorf_hombrera_pct_mediana <dbl>,
    ## #   `geomorf_espolón/gajo_pct_mediana` <dbl>,
    ## #   geomorf_vertiente_pct_mediana <dbl>,
    ## #   geomorf_vaguada_pct_mediana <dbl>,
    ## #   geomorf_piedemonte_pct_mediana <dbl>, geomorf_valle_pct_mediana <dbl>,
    ## #   geomorf_sima_pct_mediana <dbl>, Al_mediana <dbl>, B_mediana <dbl>,
    ## #   Ca_mediana <dbl>, Cu_mediana <dbl>, Fe_mediana <dbl>, K_mediana <dbl>,
    ## #   Mg_mediana <dbl>, Mn_mediana <dbl>, P_mediana <dbl>, Zn_mediana <dbl>,
    ## #   N_mediana <dbl>, N.min._mediana <dbl>, pH_mediana <dbl>,
    ## #   elevacion_media_mediana <dbl>, pendiente_media_mediana <dbl>,
    ## #   orientacion_media_mediana <dbl>, curvatura_perfil_media_mediana <dbl>,
    ## #   curvatura_tangencial_media_mediana <dbl>,
    ## #   abundancia_global_mediana <dbl>, riqueza_global_mediana <dbl>,
    ## #   heterogeneidad_ambiental_varianza <dbl>,
    ## #   geomorf_llanura_pct_varianza <dbl>, geomorf_pico_pct_varianza <dbl>,
    ## #   geomorf_interfluvio_pct_varianza <dbl>,
    ## #   geomorf_hombrera_pct_varianza <dbl>,
    ## #   `geomorf_espolón/gajo_pct_varianza` <dbl>,
    ## #   geomorf_vertiente_pct_varianza <dbl>,
    ## #   geomorf_vaguada_pct_varianza <dbl>,
    ## #   geomorf_piedemonte_pct_varianza <dbl>,
    ## #   geomorf_valle_pct_varianza <dbl>, geomorf_sima_pct_varianza <dbl>,
    ## #   Al_varianza <dbl>, B_varianza <dbl>, Ca_varianza <dbl>,
    ## #   Cu_varianza <dbl>, Fe_varianza <dbl>, K_varianza <dbl>,
    ## #   Mg_varianza <dbl>, Mn_varianza <dbl>, P_varianza <dbl>,
    ## #   Zn_varianza <dbl>, N_varianza <dbl>, N.min._varianza <dbl>,
    ## #   pH_varianza <dbl>, elevacion_media_varianza <dbl>,
    ## #   pendiente_media_varianza <dbl>, orientacion_media_varianza <dbl>,
    ## #   curvatura_perfil_media_varianza <dbl>,
    ## #   curvatura_tangencial_media_varianza <dbl>,
    ## #   abundancia_global_varianza <dbl>, riqueza_global_varianza <dbl>,
    ## #   heterogeneidad_ambiental_minimo <dbl>,
    ## #   geomorf_llanura_pct_minimo <dbl>, geomorf_pico_pct_minimo <dbl>,
    ## #   geomorf_interfluvio_pct_minimo <dbl>,
    ## #   geomorf_hombrera_pct_minimo <dbl>,
    ## #   `geomorf_espolón/gajo_pct_minimo` <dbl>,
    ## #   geomorf_vertiente_pct_minimo <dbl>, geomorf_vaguada_pct_minimo <dbl>,
    ## #   geomorf_piedemonte_pct_minimo <dbl>, geomorf_valle_pct_minimo <dbl>, …

Ejecuto también un ANOVA de una vía, de la `riqueza_global` respecto de
`habitat` de tipo `Old*` (e.g. `OldHigh`, `OldLow`, `OldSlope`)

``` r
agrupado_por_habitat %>%
  filter(str_detect(habitat, 'Old*')) %>%
  oneway.test(formula = riqueza_global ~ habitat)
```

    ## 
    ##  One-way analysis of means (not assuming equal variances)
    ## 
    ## data:  riqueza_global and habitat
    ## F = 12.572, num df = 2.000, denom df = 16.122, p-value = 0.0005125

El resultado sugiere que “existen ‘diferencias significativas’ de
`riqueza_global` entre `habitat` de tipo `Old*`”. Esta expresión,
“diferencias significativas”, resuena mucho en ciencia hoy en día, y
más bien “chirría”. De esto hablaré más adelante, por lo pronto, donde
quiera que la veas, levanta una ceja.

### `tidyr`

Te ayudará a transformar tus datos para organizarlos de forma “tidy”.
Las funciones de `tidyr` que conocerás son (también cuenta con muchas
otras): `gather()` y `spread()`.

### `readr`

Te servirá para leer eficiente y rápidamente tablas alojadas en archivos
de texto.

### `ggplot2`

Te ayudará en la visualización de tus datos, utilizando gramática de
gráficos.
