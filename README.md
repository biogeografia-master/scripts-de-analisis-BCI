# Scripts de análisis de BCI
Por José Ramón Martínez Batlle

Dentro de la asignatura ["Biogeografía (GEO-131)"](https://github.com/biogeografia-master/material-de-apoyo), de la Universidad Autónoma de Santo Domingo (UASD), asigno la redacción de un manuscrito basado en ecología numérica, utilizando como referencia el trabajo de Borcard et al. (2018).

Para obtener resultados empíricos consistentes, se necesitan datos de comunidad y ambientales de calidad. El [conjunto de datos de BCI](http://ctfs.si.edu/webatlas/datasets/bci) es idóneo en este sentido, puesto que dispone de múltiples variables ambientales e información de la comunidad basada en repetidos censos desde la década de los 80 del siglo XX (Condit, 1998; Hubbel et al., 1999; Hubbel et al., 2005). Los datos están bien documentados y contienen todas las variables necesarias para cualquier análisis de ecología numérica.

> Aviso para estudiantes. No clones este repo dentro del de manuscrito. Es preferible visualizarlo desde GitHub. Si lo clonas, asegúrate de hacerlo fuera del repo de manuscrito.

## Recursos: vídeos y scripts

### Herramientas Rmarkdown, GitHub

* [Cómo introducir referencias bibliográficas citas APA usando Bibtex, insertar imágenes y tablas externas (no generadas por R en consola) en el manuscrito formato RMarkdown](https://drive.google.com/file/d/1DK7QNUbicjoeVUusm6Kwpk299H8oD6Jj/view?usp=sharing)

* [Video corto sobre clonar un repositorio desde GitHub a RStudio, identificarme como usuario ante GitHub, hacer commit y push](https://drive.google.com/file/d/1roecMRxQYIOmUA3--pa7Syh4hQls7Y6v/view?usp=sharing)

### BCI

* [Datos censales parcela permanente 50 ha árboles BCI, explicado por el tali](https://drive.google.com/file/d/1pHQi-9NQYGQnXKXX02eZOpTDoQarvrXl/view?usp=sharing)

### Scripts de análisis

> Te recomiendo usar las versiones `*.md`  de los scripts de análisis para ver el contenido desarrollado en cada caso, y las versiones no tejidas o scripts de `*.R` crudos para copiar código, pegarlo en R y reproducir los análisis.

* Análisis exploratorio de datos
  
    * [aed_1: Crear script de análisis para generar objetos, tales como tablas de R y gráficos (por ejemplo, gráficos de mosaico, gráficos de dispersión), e insertarlos en documento de manuscrito formato RMarkdown (análisis exploratorio de datos)](aed_1.md). [Vídeo asociado. T39:19](https://drive.google.com/file/d/1gZwltgEHpF_jdzpF-hWDZWd5QyMb3iT3/view?usp=sharing)
      
    * [aed_2: Tutorial de `tidyverse`. Este vídeo tiene finalidad didáctica. Los análisis mostrados en éste son sólo demostrativos, y no tienen por qué terminar en tu manuscrito, a menos que así lo desees. En otros vídeos si ejecutaré análisis que podrían serte útiles en tu manuscrito, aplicando lo mostrado en este tutorial](aed_2_tidyverse.md). [Vídeo asociado. T1:48:09](https://drive.google.com/file/d/1YUDPMgKsxzxXSUQAImY5SvBUQdciP19G/view?usp=sharing)
      
    * [aed_3: Cómo crear mapas de riqueza y abundancia global y de mi familia de plantas asignada. En este vídeo, verás cómo crear mapas interactivos y estáticos de abundancia y riqueza, que podrás insertar en tu manuscrito](aed_3_mapas_riqueza_abun_global_mi_familia.md). [Vídeo asociado. T35:31](https://drive.google.com/file/d/1BCNyc_z6ikS2UdADzSudfoaOtCqPfX7L/view?usp=sharing)
      
    * [aed_4: Cómo crear mapas de variables ambientales](aed_4_mapas_variables_ambientales.md). [Vídeo asociado. T27:37](https://drive.google.com/file/d/1O8wfLqeTCES8M1UPi5HGuLcTqjDoNxNm/view?usp=sharing)
      
    * [aed_5: Cómo realizar análisis y paneles de correlación entre variables ambientales](aed_5_correlaciones_variables_ambientales.md). [Vídeo asociado. T40:39](https://drive.google.com/file/d/1LqxUdjn_M8qqWoIw3JqLTYcWDPnj75zE/view?usp=sharing)
      
    * [aed_6: Mapas de variables ambientales **por lotes**](aed_6_mapas_por_lotes.md). [Vídeo asociado. T15:33](https://drive.google.com/file/d/14U2Yrk4c5FelcZFznjMc0uXnBa7XfCdY/view?usp=sharing)
      
* Medición de asociación
  
    * [ma_1: Modos de análisis Q y R. Paradoja de Orlóci](medicion_asociacion_1_modo_Q_paradoja_orloci.md). [Vídeo asociado. T29:12](https://drive.google.com/file/d/1Ie3wws6ZCaKnDQtQEVMRa6TwjSni6v4J/view?usp=sharing)
      
    * [ma_2: Modo de análisis Q, para comparar asociación entre objetos (sitios) usando métricas de distancia](medicion_asociacion_2_modo_Q_mi_familia.md). [Vídeo asociado. T36:11](https://drive.google.com/file/d/12_Rha9iQFD5MWOFmwiNfXTA767d1JT2e/view?usp=sharing)
      
    * [ma_3: Modo de análisis R, para comparar dependencia entre descriptores (variables) usando tanto métricas de distancia como de correlación](medicion_asociacion_3_modo_R_mi_familia.md). [Vídeo asociado. T28:16](https://drive.google.com/file/d/1HoR1phyXfmr8cICWG36NClv17l00nyxS/view?usp=sharing)
      
* Análisis de agrupamiento (*cluster analysis*)
  
    * [aa_1: agrupamiento jerárquico](aa_analisis_de_agrupamiento_1_jerarquico.md). [Vídeo asociado. T38:14](https://drive.google.com/file/d/1nvxyygYLlaK2lFuIbfxeApGtFN85wXFU/view?usp=sharing)
      
    * [aa_2: Interpretación y comparación de resultados](aa_analisis_de_agrupamiento_2_intepretacion_resultados.md). [Vídeo asociado. T1:11:54](https://drive.google.com/file/d/14Mh6eEGJnvfy7uw969KvJGEqF7DW33rw/view?usp=sharing)
      
    * [aa_3: Grupos (clústers), variables ambientales y mapas](aa_analisis_de_agrupamiento_3_variables_ambientales_segun_grupos.md). [Vídeo asociado. T33:49](https://drive.google.com/file/d/16MpxVlnNRnYrTT7pz-ydT5n61_Ns517P/view?usp=sharing)
      
    * [aa_4: Especies indicadoras, especies con preferencia por hábitats](aa_analisis_de_agrupamiento_4_especies_indicadoras_preferencia_por_habitat.md). [Vídeo asociado. T55:05](https://drive.google.com/file/d/1EKKR7tePhnZFTButwWWRXcbgNuRnc9j3/view?usp=sharing)
      
* Análisis de ordenación simple (no restringida) y canónica (restringida)
  
    * [to_1: Ordenación no restringida. PCA, CA y PCoA](to_tecnicas_de_ordenacion_1_ordenacion_no_restringida_PCA_CA_MCA.md). [Vídeo asociado. T2:08:05](https://drive.google.com/file/d/1WpOmbkbZ4EVQBCisVb8ggDUXQ0guDmO3/view?usp=sharing)
      
    * [to_2: Ordenación restringida o ‘canónica’. RDA, CCA](to_tecnicas_de_ordenacion_2_ordenacion_restringida_RDA_CCA.md). [Vídeo asociado. T1:11:44](https://drive.google.com/file/d/18_WRG5hSxm6iftOiqMkNKKXxp3PhqBwe/view?usp=sharing)
      
* Análisis de diversidad  alpha y beta
  
    * [di_1: Diversidad alpha](di_1_analisis_de_diversidad_diversidad_alpha.md). [Vídeo asociado. T1:45:15](https://drive.google.com/file/d/1bsNkshm7PsrqkO7qMwkLdWItCeyblztE/view?usp=sharing)
      
    * [di_2: Diversidad beta](di_2_analisis_de_diversidad_diversidad_beta.md). [Vídeo asociado. T23:29](https://drive.google.com/file/d/1bsNkshm7PsrqkO7qMwkLdWItCeyblztE/view?usp=sharing)
      
* Ecología espacial
  
    * [ee: Autocorrelación T](ee_ecologia_espacial.md). [Vídeo asociado. T1:09:52](https://drive.google.com/file/d/14S4fSKfcm3T39WTLcnW69OJV2tKQtw6p/view?usp=sharing)

### Otros recursos:

* [Drive de Google conteniendo presentaciones de diapositivas](https://drive.google.com/drive/folders/12NbRrZlw6qBtaEWH6KqiAw1XaSutxFDg?usp=sharing)

 * [Introducción a R (tutorial1), para aprender sobre objetos de R](https://geofis.shinyapps.io/tutorial1/)

## Referencias

* Borcard, D., Gillet, F., & Legendre, P. (2018). *Numerical ecology with R*. Springer.

* [Hubbell, S.P., Condit, R., and Foster, R.B. 2005. Barro Colorado Forest Census Plot Data. URL http://ctfs.si.edu/webatlas/datasets/bci](http://ctfs.si.edu/webatlas/datasets/bci/)

* Condit, R. 1998. Tropical Forest Census Plots. Springer-Verlag and R. G. Landes Company, Berlin, Germany, and Georgetown, Texas.

* Hubbell, S.P., R.B. Foster, S.T. O'Brien, K.E. Harms, R. Condit, B. Wechsler, S.J. Wright, and S. Loo de Lao. 1999. Light gap disturbances, recruitment limitation, and tree diversity in a neotropical forest. Science 283: 554-557.
