#' ---
#' title: "Análisis exploratorio de datos. Colección tidyverse"
#' author: "JR"
#' date: "18 de octubre, 2020"
#' output: github_document
#' ---

#' # ¿Qué es tidyverse?
#' 
#' Es una colección de paquetes con los que podrás importar, transformar, visualizar, modelar y presentar datos. La colección se compone de 8 paquetes, de los cuales verás sobre todo 3: `dplyr`, `tidyr` y `ggplot2`.
#' 
#' Todos estos paquetes comparten estructuras comunes. Una de las herramientas que incorpora la colección es la pipa `%>%` (**SHORTCUT: `CTRL+SHIFT+M`**), la cual importa desde el paquete `magrittr`. Usarás la pipa para construir "tuberías" de procesamiento sin necesidad de crear objetos intermedios. En una tubería, puedes interpretar la pipa como **"luego"**, y verás más adelante por qué. La función principal de la pipa (tiene muchas, pero esta es la más importante) es pasar el resultado del objeto a su izquierda como primer argumento de la función a su derecha. El siguiente ejemplo explica su uso:
#' 
#' `objeto1 %>% funcion1()` es equivalente a `funcion1(argumento1 = objeto1)`
#' 
#' > La idea del *pipe* pertenece a la tradición de sistemas tipo Unix y, en origen, su función era comunicar distintos procesos, usando la salida estándar de uno (*stdout*) como entrada estándar (*stdin*) del siguiente.
#' 
#' Su ventaja radica en que, si necesitaras continuar procesando los datos, no tendrás que anidar ni crear objetos intermedios. En el siguiente ejemplo, asigno el resultado de una cadena al objeto nombrado `resultado`:
#' 
#' `resultado <- objeto1 %>% funcion1() %>% funcion2() %>% funcion3()`
#' 
#' Puedes leer lo anterior como *"objeto1 pasa como primer argumento de funcion1, **luego** el resultado de funcion1 pasa como primer argumento de funcion2, **luego** el resultado de funcion2 pasa como primer argumento de funcion3*.
#' 
#' Para replicar esta operación sin la pipa, podrías realizarlo de, por ejemplo, dos maneras distintas:
#' 
#' * Opción 1, anidar:
#' 
#' `resultado <- funcion3(funcion2(funcion1(objeto1)))`
#' 
#' Opción 2, crear objetos intermedios:
#' 
#' `tmp1 <- funcion1(objeto1)`
#' `tmp2 <- funcion2(tmp1)`
#' `resultado <- funcion3(tmp2)`
#' 
#' Notarás que, comparada con estas dos últimas opciones, la tubería es más limpia que estas dos últimas opciones. La tubería puedes leerla de forma encandenada, a diferencia del estilo anidado y de creación de objetos intermedios, que añade una cierta complejidad de lectura para el usuario/a, sobre todo para personas sin conocimientos de programación. Precisamente por esta razón fue que decidí introducir la colección tidyverse, para así mostrarte algunas ideas que podrás aplicar a tus datos y, en principio, para facilitarte la vida ("*no me ayude' tali*"). Ahora bien, si decides programar en R más adelante, deberás aprender las capacidades de programación orientada a objetos y programación funcional de R.
#' 
#' ¡Comencemos!
#' 
#' ## Paquetes
#' 
library(tidyverse)
library(sf)
#' 
#' > `sf` te ayudará a leer el objeto `bci_env_grid` como un *simple feature*, el cual se encuentra dentro del archivo `biodata/matriz_ambiental.Rdata`. Esto extenderá las capacidades espaciales del objeto.
#' 
#' ## Cargar datos
#' 
load('biodata/matriz_ambiental.Rdata')
load('biodata/Apocynaceae-Meliaceae-Sapotaceae.Rdata')
#'  
#' ## Paquete `dplyr`
#' 
#' Te servirá para manipular datos mediante verbos. Los verbos de `dplyr` que conocerás son (hay muchos otros): `select()`, `filter()`, `arrange()`, `mutate()`, `group_by()`, `summarise()` y `join`.
#' 
#' ### Verbo `select`
#' 
#' Comúnmente, necesitas seleccionar una o varias columnas de una tabla. Para esto existe el verbo `select`. Te muestro un ejemplo aplicado a la matriz de comunidad (que por ahora la verás como `simple feature`), seleccionando las columnas `id` (número identificador de quadrat) y `pH` (pH del suelo):
bci_env_grid %>%
  select(id, pH)
#' > Importante: el objeto `bci_env_grid` permanece intacto, a menos que se use dicho nombre para reasignarlo a otro objeto. Mientras no se use el asignador `<-`, sólo verás que manipulo y visualizo copias del objeto original.
#' Fíjate en la clase del objeto `bci_env_grid`. Para ello usaré la función de R `class`. No sólo se admiten verbos `dplyr`, cualquier función puede usarse:
bci_env_grid %>%
  class
#' El objeto `bci_env_grid` es a la vez de clase `sf` (*simple feature*) y `data.frame`, es decir, es tanto tabla como objeto espacial, por lo que se puede representar en un mapa. Este objeto no pierde la clase `sf`, por lo que verás que aparece información geométrica y geoespacial en el encabezado, y luego un extracto de la tabla de datos (como máximo, las 10 primeras filas). Para convertirlo a un simple `data.frame`, hay que "tumbar" su geometría con `st_drop_geometry`:
bci_env_grid %>%
  select(id, pH) %>%
  st_drop_geometry
#' Fíjate ahora en la clase de `bci_env_grid %>% select(id, pH) %>% st_drop_geometry`, que en este caso es sólo `data.frame`:
bci_env_grid %>%
  select(id, pH) %>%
  st_drop_geometry %>%
  class
#' > Al introducir un `<enter>` después de la pipa, el código puede continuar en la línea siguiente. Esto se hace para evitar que la línea de código sea legible sin necesidad de desplazarse hacia la derecha (es aconsejable no superar los 80 caracteres en una misma línea de código, según recomiendan en la [Google’s R Style Guide](https://google.github.io/styleguide/Rguide.html)). Como convención, escribiré un `<enter>` después de cada operador pipa.
#' 
#' Seleccionaré, y a la vez renombraré, dos columnas con `select` (recuerda: no estoy modificando el objeto original, simplemente trabajo en copias no asignadas). De paso, sólo mostraré las 6 primeras filas al aplicar `head` al final de la tubería (no sólo se admiten verbos `dplyr`, cualquier función de R puede entrar en la tubería):
bci_env_grid %>%
  select(id_de_quadrat = id, pH_del_suelo = pH) %>%
  st_drop_geometry %>%
  head
#' 
#' Otra funcionalidad de `select` es poder seleccionar columnas según patrones. Por ejemplo, si quisera seleccionar únicamente las variables sobre geomorfología (todas las columnas que comienzan por `geomorf`), podría hacerlo con relativa facilidad usando la función de ayuda `contains`
#' 
bci_env_grid %>%
  select(contains('geomorf')) %>%
  st_drop_geometry
#' ...y también usando expresiones regulares con `matches`, usando por ejemplo dos cadenas de caracteres...
bci_env_grid %>%
  select(matches('geomorf|habit', ignore.case = F)) %>%
  st_drop_geometry
#' ...o pidiendo todas las columnas que comienzan por mayúsculas, excepto las que comienzan por "U", lo cual excluye las coordenadas (e.g. `UTM.NS`), y deja sólo los elementos de análisis de suelo.
bci_env_grid %>%
  select(matches('^[A-T,Z]', ignore.case = F)) %>%
  st_drop_geometry
#' ### Verbo `filter`
#' 
#' Ahora mostraré sólo los elementos con `pH` mayor que 5, usando el verbo `filter`
bci_env_grid %>%
  select(id, pH) %>%
  st_drop_geometry %>%
  filter(pH>5)
#' O filtro por aquellos con `id` 31 y 50:
bci_env_grid %>%
  select(id, pH) %>%
  st_drop_geometry %>%
  filter(id == c(31, 50))
#' 
#' ### Verbo `arrange`
#' 
#' Pruebo también con la matriz de comunidad. Por ejemplo, introduzco en la tubería la función `colSums`, que devuelve un vector cuyos elementos están nombrados (tienen un atributo, en este caso, el nombre de especie), donde cada elemento representa la abundancia por especie.
mc_apcyn_melic_saptc %>%
  colSums
#' Y también obtengo la abundancia por quadrat.
mc_apcyn_melic_saptc %>%
  rowSums
#' Uso a continuación el verbo `arrange` para mostrar los registros de la matriz ambiental ordenados ascendentemente por pH.
bci_env_grid %>%
  select(id, pH) %>%
  st_drop_geometry %>%
  arrange(pH)
#' Ahora usaré `arrange` para mostrar los registros de la matriz ambiental ordenados DESCendentemente por pH.
bci_env_grid %>%
  select(id, pH) %>%
  st_drop_geometry %>%
  arrange(desc(pH))
#' 
#' ### Verbo `mutate`
#' 
#' Usaré el verbo `mutate` para crear una nueva columna. Por ejemplo, creo una columna que contenga `habitat` y `quebrada` separadas por una coma:
bci_env_grid %>%
  st_drop_geometry %>%
  select(habitat, quebrada) %>% 
  mutate(habitat_quebrada = paste(habitat, quebrada, sep = ', '))
#' Ahora `mutate`, pero con números: creo una columna de área de cada cuadro (necesitas también la función `st_area`, del paquete `sf`):
bci_env_grid %>%
  mutate(area = st_area(geometry)) %>%
  select(id, area) %>%
  st_drop_geometry %>%
  head
#' ...y ahora más complejo: obtengo la densidad de individuos por metro cuadrado, ordenados descendentemente por dicha densidad, y conservando sólo los 6 registros con mayores densidades.
bci_env_grid %>%
  mutate(area = st_area(geometry), densidad_indiv = abundancia_global/area) %>% 
  select(id, densidad_indiv) %>% 
  st_drop_geometry %>% 
  arrange(desc(densidad_indiv)) %>% 
  head
#' 
#' ### Verbos `group_by` y `summarise`
#' 
#' Los verbos `group_by` y `summarise` son útiles para producir resúmenes por grupos.
#' 
#' Agruparé la matriz ambiental por la columna `habitat`, dejando sólo las variables numericas que hagan sentido (por ejemplo, excluyo `id`, `UTM.EW`, `UTM.NS`), y lo asignaré a `agrupado_por_habitat` para luego reutilizarlo:
agrupado_por_habitat <- bci_env_grid %>%
  st_drop_geometry %>%
  group_by(habitat) %>% 
  select_if(is.numeric) %>%
  select(-id, -UTM.EW, -UTM.NS)
agrupado_por_habitat
#' Observa el encabezado: el objeto es `A tibble: 50 x 32` y hay 5 grupos (`Groups:   habitat [5]`). Calculo cuántos elementos (filas) hay por grupo:
agrupado_por_habitat %>% summarise(n = n())
#' ...y también algunos estadísticos de las columnas `pH`, `abundancia_global` y `riqueza_global` por ejemplo:
agrupado_por_habitat %>%
  summarise(
    n = n(),
    media_pH = mean(pH),
    media_abundancia = mean(abundancia_global),
    media_riqueza = mean(riqueza_global)
  )
#' ...o la media de todas las variables numéricas
agrupado_por_habitat %>%
  summarise_all(mean)
#' ...no caben, mejor por partes
agrupado_por_habitat %>%
  summarise_all(mean) %>% 
  select(1:6) %>% 
  print(width=300)
agrupado_por_habitat %>%
  summarise_all(mean) %>% 
  select(1,7:12) %>% 
  print(width=300)
agrupado_por_habitat %>%
  summarise_all(mean) %>% 
  select(1, 13:25) %>% 
  print(width=300)
agrupado_por_habitat %>%
  summarise_all(mean) %>% 
  select(1, 26:32) %>%
  print(width=300)
#' ...y no sólo un estadístico, sino varios:
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
#' Ejecuto también un ANOVA de una vía, de la `riqueza_global` respecto de `habitat` de tipo `Old*` (e.g. `OldHigh`, `OldLow`, `OldSlope`)
agrupado_por_habitat %>%
  filter(str_detect(habitat, 'Old*')) %>%
  oneway.test(formula = riqueza_global ~ habitat)
#' El resultado sugiere que "existen 'diferencias significativas' de `riqueza_global` entre `habitat` de tipo `Old*`". Esta expresión, "diferencias significativas", resuena mucho en ciencia hoy en día, y más bien "chirría". De ella hablaré más adelante, por lo pronto, donde quiera que la veas, levanta una ceja.
#' 
#' Finalmente, te muestro `join`. Más que una función, `join` es una función genérica con varios métodos para unir tablas que comparten al menos un atributo en común, siendo `inner_join` el más usado. Imagina un caso de uso: calculas abundancia de tu familia de plantas por cada quadrat de 1 ha. Ahora necesitas unir dicha información a la matriz ambiental. Aunque hay múltiples soluciones, `inner_join` es la que con mayor consistencia resuelve este problema.
#' 
#' Obtendré una tabla con dos columnas: código identificador de quadrat de 1 ha (le llamaré `id`), y abundancia de todas las plantas de mi familia por quadrat (le llamaré `abundancia_mi_familia`)
id_abundancia_fam <- mc_apcyn_melic_saptc %>%
  mutate(abundancia_mi_familia = rowSums(.)) %>% 
  rownames_to_column(var = 'id') %>%
  mutate(id = as.numeric(id)) %>% #Numérico, garantiza compatibilidad con id de bci_env_grid
  select(id, abundancia_mi_familia)
id_abundancia_fam %>% tibble
#' Dado que `id_abundancia_fam` y `bci_env_grid` comparten el campo `id`, a través de éste se puede realizar la unión. Primero escribiré `bci_env_grid`, luego la función `inner_join`, que tiene como argumentos la tabla `x` y la tabla `y`. La tabla `x`, primer argumento, será `bci_env_grid`, y la tabla `y` será `id_abundancia_fam`. El campo de unión será `id`, el cual es compartido.
bci_env_grid %>%
  inner_join(y = id_abundancia_fam, by = 'id')
#' El resultado muestra la `bci_env_grid`, ahora con los datos de mi familia como parte de la matriz. Nótese que no asigné el objeto, por lo que la matriz ambiental original sigue intacta.
#' 
#' ## `tidyr`
#' 
#' Te ayudará a reordenar datos, mediante transformación de su estructura, para organizarlos de forma "tidy". Los verbos de `tidyr` que conocerás son (también cuenta con muchas otras): `pivot_longer()` (antiguamente `gather()`) y `pivot_wider()` (antiguamente `spread()`).
#' 
#' ### Verbo `pivot_longer`
#' 
#' Cuando necesitas reunir varias columnas, o lo que es lo mismo, hacerlas que pivoten a lo largo en la tabla, utilizas este verbo. Pasas de tener una tabla organizada a lo ancho a tenerla organizada a lo largo.
#' 
#' ![](https://uc-r.github.io/public/images/dataWrangling/gather1.png)
#' *Tomado de: UC Business Analytics R Programming Guide. Reshaping Your Data with tidyr. https://uc-r.github.io/tidyr*
#' 
#' Es común realizar "reunión" columnas cuando nos interesa aplicar análisis masivos a múltiples variables utilizando un criterio de agrupamiento. 
#' 
#' Pongo un ejemplo. Por tipo de hábitat, ¿cuánto es el promedio de los porcentajes de cada uno de los 10 tipos de unidades geomorfológicas? `pivot_longer` lo hará al final, pero primero, hay que seleccionar solo las columnas de porcentajes de unidades geomorfológicas y la de habitat...
pivotpaso1 <- bci_env_grid %>%
  st_drop_geometry %>% 
  select(matches('geomorf|habitat'))
pivotpaso1 %>% tibble
#' ...luego reunir todas las columnas de geomorfología pivotando en torno a la columna `habitat`...
pivotpaso2 <- pivotpaso1 %>%
  pivot_longer(
    cols = contains('geomorf'),
    names_to = 'variable',
    values_to = 'valor')
pivotpaso2 %>% tibble
#' ...y finalmente obtener las medias de porcentajes de geomorfología por cada grupo de habitat, presentando en orden alfabético los tipos de hábitat y luego, dentro de cada tipo de hábitat, descendentemente por media.
pivotpaso3 <- pivotpaso2 %>%
  group_by(habitat, variable) %>% 
  summarise(media = mean(valor))
pivotpaso3 %>% arrange(habitat, desc(media)) %>% print(n=Inf)
#' `pivot_longer` también es útil para realizar paneles de gráficos de muchas variables, como verás en la siguiente sección.
#' 
#' La operación contraria a `pivot_longer` se realiza con `pivot_wider`. Supongamos que ahora necesitamos que la tabla anterior se muestre a lo ancho, pivtoando igualmente `habitat`, expandiendo la columna `variable` y utilizando como relleno de las celdas, los valores de la columna media. Esto se consigue así:
#' 
pivotpaso3 %>%
  ungroup() %>%
  pivot_wider(
    id_cols = habitat,
    names_from = variable,
    values_from = media)
#' 
#' ## `ggplot2`
#' 
#' Te ayudará en la visualización de tus datos, utilizando gramática de gráficos.
#'
#' Un gráfico `ggplot` utiliza capas para mostrar la información. Los objetos fuente son `data.frame`. Puedes combinar varias fuentes de datos y, de una fuente de datos, puedes estilizar cada elemento a graficar mediante las denominadas capas "estéticas" (`aes`). Igualmente, puedes combinar distintas geometrías mediante las funciones `geom_* `. La estructura `ggplot` utiliza el símbolo `+` para agregar capas.
#' 
#' Explicaré su uso con ejemplos, descomponiendo las partes de una sentencia `ggplot` para fines didácticos, aunque se pueden generar gráficos escribiendo todas las partes en una sola sentencia.
#' 
#' Primero incluiré la función `ggplot`, para crear un espacio de coordenadas según los datos disponibles en `bci_env_grid`. Lo escribo como primera parte de la estructura y, lógicamente, como no especifico datos más detalles, sólo se creará un espacio (sin coordenadas) a rellenar.
p0 <- ggplot(bci_env_grid)
p0
#' A continuación, definiré las variables estéticas sobre las que construiré la simbología, añadiéndolas al objeto anterior mediante el operador `+`. Utilizaré las variables `riqueza_global` y `abundancia_global` para hacer un gráfico de dispersión.
p1 <- p0 + aes(x = abundancia_global, y = riqueza_global)
p1
#' El espacio de coordenadas ya está creado, y `ggplot2` está preparado para aceptar geometrías. Le añadiré una capa de puntos mediante `geom_point`, a la cual no tendré que definirle coordenadas de mapeo `x` e `y`, porque las tomará de la definición global realizada en `p1`:
p2 <- p1 + geom_point()
p2
#' Dado que en `p1` definí las coordenadas de mapeo `aes(x = abundancia_global, y = riqueza_global)`, puedo añadir, además de `geom_point`, más capas que reaprovechen dicha definición. Por ejemplo, añadiré una capa de resumen a `p2` basada en modelo (`geom_smoot`), por ejemplo usando un modelo lineal `lm`.
p3 <- p2 + geom_smooth(formula = y ~ x, method = 'lm')
p3
#' En `p3`, tanto `geom_point` como `geom_smooth` aprovechan las coordenadas del mapeo definido en `p1`.
#' 
#' Una forma alterna permite definir la capa estética dentro de la geometría con resultado idéntico. Por ejemplo:
p4 <- p0 +
  geom_point(mapping = aes(x = abundancia_global, y = riqueza_global))
p4
#' Esta forma tiene la ventaja de ser más corta, pero tiene la desventaja de que impide reutilizar las coordenadas de mapeo para otras geometrías, como vimos con `geom_point` y `geom_smooth` en `p3`, que aprovechaban el mismo mapeo de coordenadas simultáneamente.
#' 
#' También definiré propiedades globales del gráfico mediante temas.
p5 <- p3 + theme_bw()
p5
p6 <- p3 + theme_classic()
p6
p7 <- p3 + theme_minimal()
p7
#' Con una variable categórica, se pueden estilizar los elementos del gráfico. Por ejemplo, haré que los puntos se coloreen en función de `habitat`.
p8 <- p0 +
  geom_point(
    mapping = aes(
      x = abundancia_global,
      y = riqueza_global,
      color = habitat))
p8
#' 
#' Ahora mostraré cómo construir el último gráfico con una sentencia de conjunto, sin reaprovechar objetos anteriores:
p9 <- ggplot(bci_env_grid) +
  geom_point(
    mapping = aes(
      x = abundancia_global,
      y = riqueza_global,
      color = habitat))
p9
#' Las posibilidades de personalización de gráficos de `ggplot2` son enormes y superan el cometido de esta introducción. Sin embargo, te mostraré algunas funcionalidades adicionales que usarás en la asignatura. Un gráfico muy usado en estadística es el diagrama de cajas (*box-plot*). Este gráfico resume, a golpe de vista, los principales estadísticos descriptivos de un conjunto de datos (cuartiles, extremos, *outliers*, dispersión). Visita [la entrada de Wikipedia sobre *box-plot](https://es.wikipedia.org/wiki/Diagrama_de_caja) para aprender más sobre este útil gráficos. Es común que el *box-plot* se use para comparar una variable entre distintos tratamientos. En el siguiente ejemplo, muestro la variable `abundancia_global` considerando a `habitat` como tratamiento:
p10 <- p0 +
  geom_boxplot(mapping = aes(x = habitat, y = abundancia_global))
p10
#' Y ejemplifico también `riqueza_global`:
p11 <- p0 +
  geom_boxplot(mapping = aes(x = habitat, y = riqueza_global))
p11
#' ...la cual muestra efectos más marcados que `abundancia_global`.
#' 
#' Los dos gráficos anteriores son muy informativos, pero tienen la desventaja de que para poder compararlos, debo moverme a lo largo de la hoja. Una solución ideal sería verlos a ambos en un panel de gráficos. Esto es posible utilizando `facet_*`, reordenando los datos previamente con `pivot_longer` de `tidyr`
#' 
#' Necesitamos tres columnas, una con los nombres de los hábitats, otra con los nombres de las variables, y otra con los valores, algo tal que...
habitat_riqueza_abundancia <- bci_env_grid %>% st_drop_geometry %>% 
  select(habitat, abundancia_global, riqueza_global) %>% 
  pivot_longer(
    cols = c(abundancia_global, riqueza_global),
    names_to = 'variable',
    values_to = 'valor')
habitat_riqueza_abundancia
#' Construiré el gráfico definiendo a `habitat` en el eje `x`, y valor en `y`, mientras que usaré `variable` para establecer facetas del panel, donde el eje `y`, que es el valor, quedará libre para cada variable, puesto que la abundancia y la riqueza tienen escalas de medición muy diferentes. Quedará algo tal que esto...
habitat_riqueza_abundancia %>% 
  ggplot() + 
  aes(x = habitat, y = valor) + 
  geom_boxplot() + 
  facet_wrap( ~ variable, scal = 'free_y')
#' En resumen, usa `tidyverse` para sacar el máximo provecho de tus datos. El paquete `dplyr` te ayuda a manipular eficientemente; `tidyr` te servirá para reordenar los datos, ya sea para obtener resúmenes estadísticos eficientemente, como para generar gráficos de resumen; `ggplot2` te ayudará con la parte gráfica, que es siempre importante en tu manuscrito.