#' ---
#' title: "Análisis exploratorio de datos. Colección tidyverse"
#' author: "JR"
#' date: "18 de octubre, 2020"
#' output: github_document
#' ---

#' ### Cargar paquetes
library(tidyverse)
library(sf)

#' ### Cargar datos
load('biodata/matriz_ambiental.Rdata')
load('biodata/Apocynaceae-Meliaceae-Sapotaceae.Rdata')

#' ### ¿Qué es tidyverse?
#' 
#' Es una colección de paquetes con los que podrás importar, transformar, visualizar, modelar y presentar datos. La colección se compone de 8 paquetes, de los cuales verás sobre todo 4: `ggplot2`, `dplyr`, `tidyr` y `readr`.
#' 
#' Todos estos paquetes comparten estructuras comunes. Una de las herramientas que incorpora la colección es la pipa `%>%` (**SHORTCUT: `CTRL+SHIFT+M`**), la cual importa desde el paquete `magrittr`. Usarás la pipa para construir "tuberías" de procesamiento sin necesidad de crear objetos intermedios. En una tubería es posible interpretar la pipa como **"luego"**, verás abajo por qué. La pipa pasa el objeto a su izquierda como primer argumento de la función a su derecha. El siguiente ejemplo explica su uso:
#' 
#' `objeto1 %>% funcion1()` es equivalente a `funcion1(argumento1 = objeto1)`
#' 
#' Su ventaja radica en que, si necesitaras continuar procesando los datos, no tendrás que anidar ni crear objetos intermedios. En el siguiente ejemplo, asigno el resultado de una cadena al objeto nombrado `resultado`:
#' 
#' `resultado <- objeto1 %>% funcion1() %>% funcion2() %>% funcion3()`
#' 
#' Puedes leer lo anterior como *"objeto1 pasa como primer argumento de funcion1, **luego** el resultado de funcion1 pasa como primer argumento de funcion2, **luego** el resultado de funcion2 pasa como primer argumento de funcion3*.
#' 
#' Para replicar esta operación sin la pipa, se necesitaría algo tal que ésto:
#' 
#' `resultado <- funcion3(funcion2(funcion1(objeto1)))`
#' 
#' O alternativamente, crear objetos intermedios:
#' 
#' `tmp1 <- funcion1(objeto1)`
#' `tmp2 <- funcion2(tmp1)`
#' `resultado <- funcion3(tmp2)`
#' 
#' Notarás que la tubería es más limpia que estas dos últimas opciones. La tubería puedes leerla de forma encandenada, mientras que el anidado y la creación de objetos intermedios añade más complejidad. Principalmente por esta razón decidí introducir la colección tidyverse para mostrarte algunas ideas que podrás aplicar a tus datos y, en principio, para facilitarte la vida (no me ayude' tali).
#' 
#' ### `dplyr`
#' 
#' Te servirá para manipular datos mediante verbos. Los verbos de `dplyr` que conocerás son (hay muchos otros): `select()`, `filter()`, `arrange()`, `mutate()`, `group_by()` y `summarise()`.
#' 
#' Algo que comúnmente haces al maniupar una tabla es seleccionar una o varias columnas. Para esto existe el verbo `select`. Te muestro un ejemplo aplicado a la matriz de comunidad, seleccionando las columnas `id` (número identificador de quadrat) y `pH` (pH del suelo):
bci_env_grid %>%
  select(id, pH)
#' > Un detalle extremadamente importante: el objeto `bci_env_grid` permanece intacto, a menos que se use dicho nombre para reasignarlo a otro objeto. Mientras no se use el asignador `<-`, sólo verás que manipulo y visualizo copias del objeto original.
#' Fíjate en la clase del objeto `bci_env_grid`:
bci_env_grid %>%
  class
#' El objeto `bci_env_grid` es a la vez de clase `sf` (*simple feature*) y `data.frame`, es decir, es tanto tabla como objeto espacial que se puede representar en un mapa. Este objeto no pierde la clase `sf`, por lo que siempre verás que aparece información geométrica y geoespacial en el encabezado, y luego un extracto de la tabla de datos (como máximo, las 10 primeras filas). Para convertirlo a un simple `data.frame`, hay que "tumbar" su geometría con `st_drop_geometry`:
bci_env_grid %>%
  select(id, pH) %>%
  st_drop_geometry
#' Fíjate ahora en la clase de `bci_env_grid %>% select(id, pH) %>% st_drop_geometry`, que en este caso es sólo `data.frame`:
bci_env_grid %>%
  select(id, pH) %>%
  st_drop_geometry %>%
  class
#' > Nota que al introducir un enter después de la pipa, el código puede continuar en la línea siguiente. Esto se hace para evitar que la línea de código sea legible sin necesidad de desplazarse hacia la derecha. Como convención, escribiré un enter después de cada operador pipa.
#' Seleccionaré y a la vez cambiaré el nombre de dos columnas con `select` (recuerda: no estoy modificando el objeto original, simplemente trabajo en copias no asignadas). De paso, sólo mostraré las 6 primeras filas al aplicar `head` al final de la tubería (no sólo se admiten verbos `dplyr`, cualquier función de R puede entrar en la tubería):
bci_env_grid %>%
  select(id_de_quadrat = id, pH_del_suelo = pH) %>%
  st_drop_geometry %>%
  head
#' 
#' Ahora mostraré sólo los elementos con pH mayor que 5, usando el verbo `filter`
bci_env_grid %>%
  select(id, pH) %>%
  st_drop_geometry %>%
  filter(pH>5)
#' O filtro por aquellos con id 31 y 50:
bci_env_grid %>%
  select(id, pH) %>%
  st_drop_geometry %>%
  filter(id == c(31, 50))
#' Pruebo también con la matriz de comunidad. Por ejemplo, introduzco en la tubería la función `colSums`, que devuelve una matriz de una fila con la abundancia por especie
mc_apcyn_melic_saptc %>%
  colSums
#' Y también obtengo la abundancia por quadrat
mc_apcyn_melic_saptc %>%
  rowSums
#' Uso a continuación el verbo `arrange` para mostrar los registros de la matriz ambiental ordenados ascendentemente por pH
bci_env_grid %>%
  select(id, pH) %>%
  st_drop_geometry %>%
  arrange(pH)
#' Ahora usaré `arrange`, para mostrar los registros de la matriz ambiental ordenados DESCendentemente por pH
bci_env_grid %>%
  select(id, pH) %>%
  st_drop_geometry %>%
  arrange(desc(pH))
#' Usaré el verbo `mutate` para crear una nueva columna. Por ejemplo, creo una columna que contenga `habitat` y `quebrada` separadas por una coma:
bci_env_grid %>%
  st_drop_geometry %>%
  select(habitat, quebrada) %>% 
  mutate(habitat_quebrada = paste(habitat, quebrada, sep = ', '))
#' Ahora algo numérico: creo una columna de área de cada cuadro (necesitas también la función `st_area`, del paquete `sf`):
bci_env_grid %>%
  mutate(area = st_area(geometry)) %>%
  select(id, area) %>%
  st_drop_geometry %>%
  head
#' ...y más complejo: obtengo la densidad de individuos por metro cuadrado, ordenados descendentemente por dicha densidad, y conservando sólo los 6 registros con mayores densidades.
bci_env_grid %>%
  mutate(area = st_area(geometry), densidad_indiv = abundancia_global/area) %>% 
  select(id, densidad_indiv) %>% 
  st_drop_geometry %>% 
  arrange(desc(densidad_indiv)) %>% 
  head
#' Finalmente, te muestro los verbos `group_by` y `summarise`, los cuales son útiles para producir resúmenes por grupos.
#' Agruparé la matriz ambiental por la columna `habitat`, dejando sólo las variables numericas que hagan sentido (por ejemplo, excluyo `id`, `UTM.EW`, `UTM.NS`), y lo asignaré a `agrupado_por_habitat`:
agrupado_por_habitat <- bci_env_grid %>%
  st_drop_geometry %>%
  group_by(habitat) %>% 
  select_if(is.numeric) %>%
  select(-id, -UTM.EW, -UTM.NS)
agrupado_por_habitat
#' Observa el encabezado: el objeto es `A tibble: 50 x 32` y hay 5 grupos (`Groups:   habitat [5]`). Obtendré cuántos elementos (filas) hay por grupo:
agrupado_por_habitat %>% summarise(n = n())
#' ...y también los estadísticos de las columnas `pH`, `abundancia_global` y `riqueza_global` por ejemplo:
agrupado_por_habitat %>%
  summarise(
    n = n(),
    media_pH = mean(pH),
    media_abundancia = mean(abundancia_global),
    media_riqueza = mean(riqueza_global)
  )
#' ...o los mismos estadísticos de todas las variables numéricas
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
#' Ejecuto también un ANOVA de una vía, de la `riqueza_global` respecto de `habitat` de tipo `Old*` (e.g. `OldHigh`, `OldLow`, `OldSlope`)
agrupado_por_habitat %>%
  filter(str_detect(habitat, 'Old*')) %>%
  oneway.test(formula = riqueza_global ~ habitat)
#' El resultado sugiere que "existen 'diferencias significativas' de `riqueza_global` entre `habitat` de tipo `Old*`". Esta expresión, "diferencias significativas", resuena mucho en ciencia hoy en día, y más bien "chirría" hoy. De esto hablaré más adelante, por lo pronto, donde quiera que la veas, levanta una ceja.

#' ### `tidyr`
#' 
#' 
#' Te ayudará a transformar tus datos para organizarlos de forma "tidy". Las funciones de `tidyr` que conocerás son (también cuenta con muchas otras): `gather()` y `spread()`.
#' 
#' ### `readr`
#' 
#' Te servirá para leer eficiente y rápidamente tablas alojadas en archivos de texto.
#' 
#' ### `ggplot2`
#' 
#' Te ayudará en la visualización de tus datos, utilizando gramática de gráficos.

