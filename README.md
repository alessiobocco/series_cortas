# Scripts en R para el cálculo de índices de sequía en estaciones automáticas (EMAs)

Este repositorio contiene código para calcular índices de sequía estandarizados (SPI, SPEI) para estaciones meteorológicas con series temporales cortas derivadas de estaciones automáticas. Se han implementado tres metodologías diferentes para este cálculo:

 * **Asignación del Vecino Más Cercano**: Si la estación vecina está dentro de un rango de distancia determinado por el usuario, se asignan los valores de la estación vecina más cercana.
 * **Ponderación por el Inverso de la Distancia**: Se ponderan los efectos de los N vecinos más cercanos utilizando el inverso de la distancia.
 * **Interpolación usando Kriging**: Se utiliza Kriging para estimar los parámetros de los índices para estaciones más alejadas. Este método permite la estimación de la variabilidad de los parámetros y generar N realizaciones de índices de sequía para estimar la variabilidad o el error de estimación.
   
Contenido del Repositorio

* data: Ubicación donde se encuentran los datos de entrada.
* EstadísticasMoviles: Contiene los scripts necesarios para la agregación diferentes escalas temporales. Además, se encuentra un documento que explica el mecanismo de agregación de variables (precipitación y temperaturas máxima y mínima) utilizando ventanas móviles basadas en péntadas.
* IndicesSequia: Contiene los scripts para el cálculo de índices de sequía en alta frecuencia. Además, se incluye un documento: a) que detalla los índices que se utilizan para el análisis de las sequías y b) caracterización de los eventos secos. Se describe para cada uno de los índices su forma de cálculo.
* IdentificarEventos: Contiene los scripts para la identificación de eventos secos o humedos a partir de indices de sequía estandarizados (SPI y SPEI). Además, se incluye un documento que  define el concepto de un evento seco, como identificarlo y las métricas que se utilizan para caracterizar cada evento.

