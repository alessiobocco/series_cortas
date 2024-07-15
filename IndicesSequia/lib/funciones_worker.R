# -----------------------------------------------------------------------------#
# --- Funciones a ejecutar por parte de cada worker ----
# -----------------------------------------------------------------------------#

CalcularIndicesSequiaUbicacion <- function(input.value, script, config, configuraciones.indices, estadisticas.moviles, estaciones, indice_parametro) {
  # Obtener la ubicación para la cual se calcularán los índices
  ubicacion <- input.value
  
  # Identificar la columna con el id de la ubicación (usualmente station_id, o point_id)
  id_column <- IdentificarIdColumn(ubicacion)
  
  # Informar estado de la ejecución
  script$info(glue::glue("Procesando estadisticas para la ubicación con ",
                         "{id_column} = {ubicacion %>% dplyr::pull(!!id_column)} ",
                         "(lon: {ubicacion$lon_dec}, lat: {ubicacion$lat_dec})"))
  
  # Asignar funcion rgenlog a Global Environment (sino la ejecucion en procesos hijos no funciona bien)
  assign("rgenlog", rgenlog, .GlobalEnv)
  
  
  #############################################################################################
  ## Definición de estadisticas.variables.completas, fecha.primeros.datos y fecha.ultimos.datos
  
  # Filtrar estadisticas.moviles para quedarnos solo con las asociadas a la ubicación analizada
  estadisticas.moviles.ubicacion = estadisticas.moviles %>% 
    dplyr::filter(!!rlang::sym(id_column) == dplyr::pull(ubicacion, !!id_column))
  
  # Estadísticas móviles para prcp
  estadisticas.precipitacion <- estadisticas.moviles.ubicacion %>% 
    dplyr::filter(variable_id == 'prcp', estadistico == 'Suma') %>%
    dplyr::select(fecha_desde, fecha_hasta, ancho_ventana_pentadas, prcp = valor) 
  # Estadísticas móviles para tmin
  estadisticas.temp.min      <- estadisticas.moviles.ubicacion %>% 
    dplyr::filter(variable_id == 'tmin', estadistico == 'Media') %>%
    dplyr::select(fecha_desde, fecha_hasta, ancho_ventana_pentadas, tmin = valor) 
  # Estadísticas móviles para tmax
  estadisticas.temp.max      <- estadisticas.moviles.ubicacion %>% 
    dplyr::filter(variable_id == 'tmax', estadistico == 'Media') %>%
    dplyr::select(fecha_desde, fecha_hasta, ancho_ventana_pentadas, tmax = valor) 
  
  rm(estadisticas.moviles.ubicacion); invisible(gc())
  
  # Calcular fecha de ultimos datos y alinear series de variables climaticas por fechas. 
  # Agregar columna con dato de evapotranspiracion potencial calculada a partir de
  # metodo de Hargreaves-Samani. Este metodo solamente requiere valores de precipitacion,
  # temperatura minima, maxima y la latitud del punto. 
  #
  # El problema es que el paquete que utiliza constantes calibradas para calculos a nivel 
  # mensual (6 pentadas). Por lo tanto, hay que hacer los calculos mensuales, interpolar el
  # valor de et0 (que es ciclico) para completar faltantes y luego agregar a nivel de mayor
  # cantidad de meses.
  estadisticas.variables <- estadisticas.precipitacion %>%
    dplyr::left_join(estadisticas.temp.min, by = c("fecha_desde", "fecha_hasta", "ancho_ventana_pentadas")) %>%
    dplyr::left_join(estadisticas.temp.max, by = c("fecha_desde", "fecha_hasta", "ancho_ventana_pentadas")) %>%
    dplyr::mutate(srad = CalcularRadiacionSolarExtraterrestre(fecha_desde, fecha_hasta, ubicacion$lat_dec)) %>%
    dplyr::mutate(et0 = SPEI::hargreaves(Tmin = tmin, Tmax = tmax, Pre = prcp, Ra = srad, na.rm = TRUE)) %>%
    dplyr::select(-srad, -tmin, -tmax)
  
  fecha.primeros.datos   <- min(estadisticas.variables$fecha_desde)
  fecha.ultimos.datos    <- max(estadisticas.variables$fecha_hasta)
  
  # Se borran datos que ya no será utilizados
  rm(estadisticas.precipitacion, estadisticas.temp.max, estadisticas.temp.min); invisible(gc())
  
  # Interpolacion de et0 a nivel mensual para completar faltantes
  estadisticas.mensuales.variables <- estadisticas.variables %>%
    dplyr::filter(ancho_ventana_pentadas == 6) %>%
    dplyr::mutate(pentada_inicio = fecha.a.pentada.ano(fecha_desde))
  et0        <- dplyr::pull(estadisticas.mensuales.variables, et0)
  pentada    <- dplyr::pull(estadisticas.mensuales.variables, pentada_inicio)
  tryCatch({
    fit.et0    <- mgcv::gam(et0 ~ s(pentada, bs = "cc"), method = "REML", na.action = na.omit)
    pent.pred  <- sort(unique(pentada))
    et0.pred   <- mgcv::predict.gam(object = fit.et0, newdata = data.frame(pentada = pent.pred))
    et0.smooth <- data.frame(pentada_inicio = pent.pred, et0_pred = et0.pred)
    estadisticas.mensuales.et0 <- estadisticas.mensuales.variables %>%
      dplyr::inner_join(et0.smooth, by = "pentada_inicio") %>%
      dplyr::mutate(et0_completo = dplyr::if_else(! is.na(et0), et0, et0_pred)) %>%
      dplyr::select(realizacion, fecha_desde, fecha_hasta, ancho_ventana_pentadas, et0_completo)
    rm(fit.et0, pent.pred, et0.pred, et0.smooth)
    
    # Ahora, elimino el valor de et0 de todos los niveles y dejo solo el nivel mensual
    estadisticas.variables <- estadisticas.variables %>%
      dplyr::left_join(estadisticas.mensuales.et0, by = c("realizacion","fecha_desde", "fecha_hasta", "ancho_ventana_pentadas")) %>%
      dplyr::select(-et0) %>%
      dplyr::rename(et0 = et0_completo)
    rm(estadisticas.mensuales.et0)
    
    # Por ultimo, para los valores faltantes de et0 cuyo ancho de ventana sea mayor a 6 pentadas,
    # sumo los valores correspondientes a los subgrupos de 6 pentadas que componen el periodo.
    estadisticas.variables.completas <<- purrr::map_dfr(
      .x = unique(estadisticas.variables$realizacion),
      .f = function(r){
        estadisticas.variables %>% dplyr::filter(realizacion == r) %>%
          dplyr::mutate(et0_completo = AgregarET0(fecha_desde, fecha_hasta, ancho_ventana_pentadas, et0)) %>%
          dplyr::select(-et0) %>% dplyr::rename(et0 = et0_completo)
      })
  }, error = function(e) {
    script$warn(glue::glue("No es posible ajustar ciclo estacional para ET0: {e$message}\n"))
    estadisticas.variables.completas <<- estadisticas.variables
  })
  rm(et0, pentada, estadisticas.mensuales.variables, estadisticas.variables); invisible(gc())
  
  
  #########################################################################
  ## Definición de fechas.procesables, fechas.pentada.ano y pentadas.unicas
  
  # Procesar desde la primera pentada con datos estadisticos.
  fechas.procesables <- seq.pentadas(fecha.primeros.datos, fecha.ultimos.datos)  %>%
    purrr::flatten_dbl(.) %>%
    as_date(., origin = lubridate::origin)
  
  # Determinar que fechas pertenecen a la misma pentada del ano.
  # Esas fechas tiene los mismos parametros. De este modo, optimizamos calculos.
  fechas.pentada.ano <- fecha.a.pentada.ano(fechas.procesables)
  pentadas.unicas    <- sort(unique(fechas.pentada.ano))
  
  
  ################################
  ## Inicio del cálculo de índices
  
  # Ejecutar calculo para cada configuracion.
  # Cada configuracion es una combinacion unica de indice, escala, distribucion,
  # metodo de ajuste y periodo de referencia.
  resultado <- purrr::map_dfr(
    .x = seq(from = 1, to = nrow(configuraciones.indices)),
    .f = function(row_index) {
      # Obtener configuracion de calculo
      configuracion.indice <- configuraciones.indices[row_index, ]
      
      # Informar estado de la ejecución
      script$info(glue::glue("Procesando estadisticas para la ubicación: ",
                             "{ubicacion %>% dplyr::pull(!!id_column)}, config: {configuracion.indice$id} ",
                             "(indice = {configuracion.indice$indice}, escala = {configuracion.indice$escala}, ",
                             "distribucion = {configuracion.indice$distribucion})"))
      
      # Obtener estadisticas para esa escala de tiempo
      estadisticas.variables <- estadisticas.variables.completas %>%
        dplyr::filter(ancho_ventana_pentadas == (configuracion.indice$escala * 6))
      
      #### ----------------------------------------------------------------------#
      #### Para cada fecha procesable:
      #### 1. Buscar parametros de ajuste para esta configuracion, red, ubicacion,
      ####    y pentada de fin (correspondiente a la fecha procesable)
      #### 2. Si no exiten parametros:
      ####    a. Determinar periodo de agregacion. El final del periodo esta dado
      ####       por la fecha procesable. El inicio, de acuerdo a la longitud de la escala
      ####    b. Agregar valores para el periodo de referencia. Cuando haya periodos que caigan
      ####       parte en un ano y en otro, siempre la extension es hacia el fin del periodo
      ####       de referencia asumiendo que es mas probable que haya datos mas recientes.
      ####    c. Ajustar distribucion segun parametros de configuracion. En este paso deben
      ####       considerarse la cantidad de faltantes. Si hay demasiados faltantes (de acuerdo
      ####       a unbrales definidos), entonces todos los parametros devueltos son NA.
      ####    d. Aplicar tests de bondad de ajuste (solo para el caso de ajustes parametricos
      ####       con distribucion definina). Si alguno de los tests falla, entonces considerar 
      ####       que el ajuste no es bueno. En ese caso, reemplazar todos los valores de parametros por NA.
      ####    e. Guardar parametros en base de datos. En el caso de ser un objeto de ajuste, 
      ####       no se guarda en la base de datos. Ir directamente al paso 3.
      #### 3. Una vez obtenidos los parametros de ajuste:
      ####    i. Si alguno de los parametros es NA, tanto el valor del indice como su
      ####       percentil asociado son NA. Ir directamente al paso iii.
      ####    ii. Calcular valor de indice y percentil asociado.
      ####    iii. Guardar valor de indice y percentil asociado en base de datos.
      #### ----------------------------------------------------------------------#
      
      #script$info(glue::glue("Calculando indices de sequia de la ubicación {ubicacion$point_id} ",
      #                       "para la configuración {configuracion.indice$id}"))
      
      resultados.indice.configuracion <- purrr::map_dfr(
        .x = pentadas.unicas,
        .f = function(pentada.ano) {
          # a. Buscar fechas procesables correspondientes a esa pentada
          fechas.procesables.pentada <- fechas.procesables[which(fechas.pentada.ano == pentada.ano)]
          
          # Identificar estaciones vecinas más cercanas y asignar
          # parametros a la ubicacion sobre la que se esta iterando.
          if (!is.na(config$params$vecino.mas.cercano$distancia)) {
            
            # Informar estado de la ejecución
            script$info(glue::glue("Identificando vecinos para la ubicación: ",
                                   "{ubicacion %>% dplyr::pull(!!id_column)} ",
                                   "en un rango de {config$params$vecino.mas.cercano$distancia} km"))
            
            # Convertir a matriz de coordenadas
            coords <- estaciones %>% dplyr::select(lon_dec, lat_dec) %>% as.matrix()
            # Calcular distancias en metros desde la ubicación objetivo a cada estación
            distancias <- geosphere::distVincentyEllipsoid(c(ubicacion$lon_dec, ubicacion$lat_dec), coords)
            # Convertir distancias a kilómetros
            distancias_km <- distancias / 1000
            # Agregar las distancias al dataframe
            estaciones_cercanas <- estaciones %>%
              dplyr::mutate(distancia_objetivo = distancias_km) %>%
              # Filtrar estaciones dentro del rango de N km
              dplyr::filter(distancia_objetivo <= config$params$vecino.mas.cercano$distancia)
            
            parametros.vecinos <- indice_parametro %>%
              dplyr::filter(indice_configuracion_id == configuracion.indice$id,
                            pentada_fin == pentada.ano,
                            omm_id %in% estaciones_cercanas$omm_id)
            
            if (nrow(parametros.vecinos) == 0) {
              # Informar estado de la ejecución
              script$info(glue::glue("No hay vecinos para la ubicación: ",
                                     "{ubicacion %>% dplyr::pull(!!id_column)} ",
                                     "en un rango de {config$params$vecino.mas.cercano$distancia} km, ",
                                     "se usarán parámetros interpolados"))
            }
            
            if (nrow(parametros.vecinos) > 0) {
              if (config$params$vecino.mas.cercano$ponderacion) {
                
                parametros.vecinos %<>%
                  dplyr::left_join(estaciones_cercanas %>% dplyr::select(omm_id, distancia_objetivo)) %>%
                  dplyr::group_by(parametro, indice_configuracion_id, pentada_fin) %>%
                  dplyr::mutate(
                    pesos_invertidos = 1 / distancia_objetivo,
                    pesos_invertidos_normalizados = pesos_invertidos / sum(pesos_invertidos)
                  ) %>%
                  dplyr::summarise(valor = weighted.mean(valor, pesos_invertidos_normalizados)) %>%
                  dplyr::ungroup() %>%
                  dplyr::mutate(station_id = ubicacion$station_id) %>%
                  dplyr::select(indice_configuracion_id, station_id, pentada_fin, parametro, valor)
                
              } else {
                
                parametros.vecinos %<>%
                  dplyr::left_join(estaciones_cercanas %>% dplyr::select(omm_id, distancia_objetivo)) %>%
                  dplyr::group_by(parametro, indice_configuracion_id, pentada_fin) %>%
                  dplyr::arrange(distancia_objetivo) %>%
                  dplyr::slice(1) %>%
                  dplyr::ungroup() %>%
                  dplyr::mutate(station_id = ubicacion$station_id) %>%
                  dplyr::select(indice_configuracion_id, station_id, pentada_fin, parametro, valor)
              }
            }
          }
          
          # Usar parametros de vecinos si existen
          if (nrow(parametros.vecinos) > 0 & !config$params$simulacion.parametros$simular) {
            parametros.ajuste <- parametros.vecinos %>%
              dplyr::mutate(realizacion = 1)
          } else {
            
            if (config$params$simulacion.parametros$simular) {
              
              parametros.ajuste <- purrr::map_dfr(
                .x = 1:config$params$simulacion.parametros$n_muestras,
                .f = function(realizacion) {
                  feather::read_feather(glue::glue("{config$dir$data}/{config$files$indices_sequia$parametros_intermedios}_{ubicacion$station_id}_{configuracion.indice$internal_id}_{pentada.ano}_{realizacion}.feather")) %>%
                    dplyr::mutate(realizacion)
                }
              )
            } else {
              parametros.ajuste <- feather::read_feather(glue::glue("{config$dir$data}/{config$files$indices_sequia$parametros_intermedios}_{ubicacion$station_id}_{configuracion.indice$internal_id}_{pentada.ano}.feather")) %>%
                dplyr::mutate(realizacion = 1)
            }
          }
          
          # c. Calcular indices para fechas correspondientes a esa pentada (los parametros son los mismos)
          resultados.indice.configuracion.pentada <- purrr::map_dfr(
            .x = fechas.procesables.pentada,
            .f = function(fecha.procesable) {
              purrr::map_dfr(
                .x = unique(parametros.ajuste$realizacion),
                .f = function(r) {
                  parametros.ajuste.realizacion <- parametros.ajuste %>% dplyr::filter(realizacion == r)
                  return (tibble::tibble(realizacion = r) %>% tidyr::crossing(
                    CalcularIndicesSequiaUbicacionFecha(ubicacion, fecha.procesable, parametros.ajuste = parametros.ajuste.realizacion, 
                                                        configuracion.indice, 
                                                        script, config,estadisticas.variables)))
                }
              )
            }
          )
          return (resultados.indice.configuracion.pentada)
        }
      )
      
      # iii. Guardar valor de indice y percentil asociado en base de datos.
      if (nrow(resultados.indice.configuracion) > 0) {
        script$info(glue::glue("Retornando indices de sequia de la ubicación {ubicacion %>% dplyr::pull(!!id_column)} ",
                               "para la configuración {configuracion.indice$id}"))
        valores.indice <- resultados.indice.configuracion %>%
          dplyr::mutate(conf_id = configuracion.indice$id, indice = configuracion.indice$indice, 
                        escala = configuracion.indice$escala, distribucion = configuracion.indice$distribucion, 
                        metodo_ajuste = configuracion.indice$metodo_ajuste, 
                        referencia_comienzo = as.Date(configuracion.indice$referencia_comienzo), 
                        referencia_fin = as.Date(configuracion.indice$referencia_fin),
                        metodo_imputacion_id = 0, !!id_column := dplyr::pull(ubicacion, !!id_column)) %>%
          dplyr::select(!!id_column, realizacion, pentada_fin, ano, metodo_imputacion_id, conf_id, indice, escala, distribucion,
                        metodo_ajuste, referencia_comienzo, referencia_fin, valor_dato, valor_indice, percentil_dato)
        return(valores.indice)
      } else {
        type_of_id_col <- typeof(dplyr::pull(ubicacion,!!id_column))
        valores.indice <- tibble::tibble(!!id_column := if(type_of_id_col == "integer") integer() else 
          if(type_of_id_col == "numeric") double() else 
            if(type_of_id_col == "logical") logical() else 
              if(type_of_id_col == "character") character() else
                character(), 
          realizacion = integer(), pentada_fin = double(), ano = double(), metodo_imputacion_id = double(), 
          conf_id = integer(), 
          indice = character(), escala = integer(), distribucion = character(), 
          metodo_ajuste = character(), referencia_comienzo = as.Date(character()), 
          referencia_fin = as.Date(character()),
          valor_dato = double(), valor_indice = double(), percentil_dato = double())
        return(valores.indice)
      }
    }
  )
  
  return (resultado)
}

AjustarParametrosUbicacionFecha <- function(ubicacion, fecha.procesable, configuracion.indice,
                                            script, config, estadisticas.variables, id_column) {
  # 1. Buscar parametros de ajuste para esta configuracion, ubicacion,
  #    y pentada de fin (correspondiente a la fecha procesable)
  ano.fin     <- lubridate::year(fecha.procesable)
  pentada.fin <- fecha.a.pentada.ano(fecha.procesable)
  
  script$info(glue::glue("Ajustando parametros para pentada {pentada.fin} de la ubicación ", 
                         "{ubicacion %>% dplyr::pull(!!id_column)} y la configuración {configuracion.indice$id}"))
  
  # a. Buscar valores para el periodo de referencia. Cuando el indice es SPEI,
  #    el valor buscado es (prcp - et0), sino es (prcp)
  anos.periodo.referencia <- seq(from = lubridate::year(configuracion.indice$referencia_comienzo),
                                 to = lubridate::year(configuracion.indice$referencia_fin))
  fechas.fin.pentada.fin  <- fecha.fin.pentada(pentada.ano.a.fecha.inicio(pentada.ano = pentada.fin, ano = anos.periodo.referencia))
  datos.referencia        <- estadisticas.variables %>%
    dplyr::filter(fecha_hasta %in% fechas.fin.pentada.fin)
  if (configuracion.indice$indice == "SPEI") {
    variable.acumulada <- datos.referencia %>%
      dplyr::mutate(diferencia = prcp - et0) %>%
      dplyr::filter(! is.na(diferencia)) %>% 
      dplyr::pull(diferencia)
  } else {
    variable.acumulada <- datos.referencia %>%
      dplyr::filter(! is.na(prcp)) %>% 
      dplyr::pull(prcp)
  }
  rm(datos.referencia)
  
  # c. Ajustar distribucion segun parametros de configuracion. En este paso deben
  #    considerarse la cantidad de faltantes. Si hay demasiados faltantes (de acuerdo
  #    a unbrales definidos), entonces todos los parametros devueltos son NA.
  min.valores.necesarios <- round(length(anos.periodo.referencia) * config$params$min.proporcion.disponibles.referencia)
  configuracion.ajuste   <- config$params$ajuste.general[[configuracion.indice$metodo_ajuste]]
  funcion.ajuste         <- configuracion.ajuste$funcion
  argumentos.ajuste      <- configuracion.ajuste$parametros
  argumentos.indice      <- config$params$ajuste.particular[[configuracion.indice$metodo_ajuste]]
  if (! is.null(argumentos.indice)) {
    argumentos.indice    <- argumentos.indice[[configuracion.indice$indice]]
    if (! is.null(argumentos.indice)) {
      argumentos.ajuste  <- append(argumentos.ajuste, argumentos.indice)
    }
  }
  if ("distribucion" %in% names(formals(funcion.ajuste))) {
    argumentos.ajuste$distribucion <- configuracion.indice$distribucion
  }
  if (length(variable.acumulada) >= min.valores.necesarios) {
    argumentos.ajuste$x  <- variable.acumulada
  } else {
    script$warn(paste("Verificar los resultados, los parámetros parecen no haberse generado correctamente!!",
                      "No se pudo completar 30 observaciones, probablemente la combinación de años y",
                      "realizaciones no llega a 30 observaciones!!"))
  }
  parametros.ajuste      <- do.call(what = funcion.ajuste, args = argumentos.ajuste)
  
  # d. Aplicar tests de bondad de ajuste (solo para el caso de ajustes parametricos
  #    con distribucion definina). Si alguno de los tests falla, entonces considerar 
  #    que el ajuste no es bueno. En ese caso, reemplazar todos los valores de parametros por NA.
  if (! is.na(configuracion.indice$distribucion) || 
      (configuracion.indice$metodo_ajuste %in% c("NoParametrico", "Empirico"))) {
    script$info(glue::glue("Verificando bondad de ajuste para pentada {pentada.fin} de la ubicación ", 
                           "{ubicacion %>% dplyr::pull(!!id_column)} y la configuración {configuracion.indice$id}"))
    
    parametros.test <- list(x = variable.acumulada, config = config)
    if (configuracion.indice$metodo_ajuste == "NoParametrico") {
      parametros.test$objeto.ajuste <- parametros.ajuste  
    } else if (configuracion.indice$metodo_ajuste == "Empirico") {
      parametros.test$percentiles.ajuste <- unlist(parametros.ajuste)
    } else {
      parametros.test$distribucion      <- configuracion.indice$distribucion
      parametros.test$parametros.ajuste <- parametros.ajuste
    }
    resultado.tests <- do.call(what = "TestearBondadAjuste", args = parametros.test)
    estadisticos    <- resultado.tests$estadisticos
    if (! resultado.tests$pasa.tests) {
      # Agregar los parametros originales a la tabla de resultados de tests.
      # Solamente aplicable para ajuste parametrico
      if (configuracion.indice$metodo_ajuste != "NoParametrico") {
        estadisticos <- rbind(estadisticos, ParametrosADataFrame(parametros.ajuste) %>%
                                dplyr::mutate(test = "") %>%
                                dplyr::select(test, parametro, valor)
        )
      }
      
      # Luego, poner todos los parametros originales en NA
      parametros.na        <- as.list(rep(NA, length(parametros.ajuste)))
      names(parametros.na) <- names(parametros.ajuste)
      parametros.ajuste    <- parametros.na
    }
    
    # Guardar resultados de tests en base de datos
    if (! is.null(estadisticos)) {
      resultado.tests <- estadisticos %>%
        dplyr::mutate(!!id_column := dplyr::pull(ubicacion, !!id_column), 
                      indice_configuracion_id = configuracion.indice$id, 
                      metodo_imputacion_id = 0, pentada_fin = pentada.fin) %>%  
        dplyr::select(indice_configuracion_id, !!id_column, pentada_fin, test, parametro, metodo_imputacion_id, valor)
      feather::write_feather(resultado.tests, glue::glue("{config$dir$data}/control/resultados_tests/",
                                                         "resultados_tests_{ubicacion %>% dplyr::pull(!!id_column)}_",
                                                         "{configuracion.indice$id}_{pentada.fin}.feather"))
    }
  }
  
  # e. Guardar parametros en base de datos. En el caso de ser un objeto de ajuste, 
  #    no se guarda en la base de datos. Ir directamente al paso 3.
  if (is.list(parametros.ajuste) && (class(parametros.ajuste) != "logspline")) {
    parametros.ajuste <- ParametrosADataFrame(parametros.ajuste) %>%
      dplyr::mutate(!!id_column := dplyr::pull(ubicacion, !!id_column),  
                    indice_configuracion_id = configuracion.indice$id, 
                    metodo_imputacion_id = 0, pentada_fin = pentada.fin) %>%
      dplyr::select(indice_configuracion_id, !!id_column, pentada_fin, parametro, metodo_imputacion_id, valor)
    feather::write_feather(parametros.ajuste, glue::glue("{config$dir$data}/control/parametros/",
                                                         "parametros_{ubicacion %>% dplyr::pull(!!id_column)}_",
                                                         "{configuracion.indice$id}_{pentada.fin}.feather"))
  }
  
  return (parametros.ajuste)
}


CalcularIndicesSequiaUbicacionFecha <- function(ubicacion, fecha.procesable, parametros.ajuste, configuracion.indice, 
                                                script, config, estadisticas.variables) {
  # Calculo de pentadas de inicio y fin
  # 1. Buscar parametros de ajuste para esta configuracion, ubicacion,
  #    y pentada de fin (correspondiente a la fecha procesable)
  ano.fin               <- lubridate::year(fecha.procesable)
  pentada.fin           <- fecha.a.pentada.ano(fecha.procesable)
  fecha.fin.pentada.fin <- fecha.fin.pentada(fecha.procesable)
  
  # 2. Una vez obtenidos los parametros de ajuste:
  #    i. Si alguno de los parametros es NA, tanto el valor del indice como su
  #       percentil asociado son NA. Ir directamente al paso iii.
  calcular.indice <- TRUE
  if (is.data.frame(parametros.ajuste)) {
    # Si no hay ningun parametro nulo, pasar los parametros a lista.
    calcular.indice   <- ! any(is.na(parametros.ajuste$valor))
    parametros.ajuste <- ParametrosALista(parametros.ajuste)
  } else {
    # Es el objeto de ajuste devuelto por "logspline"
    calcular.indice <- (class(parametros.ajuste) == "logspline")
  }
  
  # ii. Calcular valor de indice y percentil asociado.
  #     Buscar primero el valor de correspondiente a la fecha procesable.
  estadisticas.calculo  <- estadisticas.variables %>%
    dplyr::filter(fecha_hasta == fecha.fin.pentada.fin)
  var.acumulada.calculo <- NA  
  if (configuracion.indice$indice == "SPEI") {
    if (nrow(estadisticas.calculo) > 0) {
      # Restarle la et0
      var.acumulada.calculo <- estadisticas.calculo %>%
        dplyr::mutate(diferencia = prcp - et0) %>%
        dplyr::pull(diferencia)
    }
  } else {
    if (nrow(estadisticas.calculo) > 0) {
      var.acumulada.calculo <- estadisticas.calculo %>%
        dplyr::pull(prcp)
    }
  }
  rm(estadisticas.calculo)
  
  # Ejecutar calculo si los parametros son validos y el valor de entrada no es NA
  if (calcular.indice && ! is.na(var.acumulada.calculo)) {
    configuracion.calculo <- config$params$calculo[[configuracion.indice$indice]]
    funcion.calculo       <- configuracion.calculo$funcion
    parametros.calculo    <- configuracion.calculo$parametros.adicionales
    parametros.calculo$x  <- var.acumulada.calculo
    if (class(parametros.ajuste) == "logspline") {
      if (! is.null(configuracion.calculo$parametro.logspline)) {
        parametros.calculo[[configuracion.calculo$parametro.logspline]] <- parametros.ajuste
      }
    } else if (is.list(parametros.ajuste)) {
      if (! is.null(configuracion.calculo$parametro.lista)) {
        parametros.calculo[[configuracion.calculo$parametro.lista]] <- parametros.ajuste
      } else if (! is.null(configuracion.calculo$parametro.vector)) {
        parametros.calculo[[configuracion.calculo$parametro.vector]] <- unname(unlist(parametros.ajuste))
      }
    }
    
    resultado.calculo <- do.call(what = funcion.calculo, args = parametros.calculo)
    return (data.frame(ano = ano.fin, pentada_fin = pentada.fin, valor_dato = resultado.calculo$valor_dato, 
                       valor_indice = resultado.calculo$valor_indice, percentil_dato = resultado.calculo$percentil_dato))
  } else {
    return (data.frame(ano = ano.fin, pentada_fin = pentada.fin, valor_dato = var.acumulada.calculo, valor_indice = NA, percentil_dato = NA))
  }
}

# Función para estimar parámetros para una combinación dada
EstimarParametrosIndices <- function(indice, pentada, param, estaciones, indice_parametro, grilla_regular, config) {
  indice_parametro_i <- indice_parametro %>%
    dplyr::filter(indice_configuracion_id == indice,
                  pentada_fin == pentada, 
                  parametro == param) %>%
    dplyr::left_join(estaciones, by = "omm_id") %>%
    tidyr::drop_na() %>%
    dplyr::select(omm_id, lon_dec, lat_dec, pentada_fin, valor) %>%
    sf::st_as_sf(coords = c("lon_dec", "lat_dec"), crs = config$proj4string$latlon) %>%
    sf::st_transform(crs = config$proj4string$planar) %>%
    dplyr::mutate(
      x = sf::st_coordinates(.)[, "X"],
      y = sf::st_coordinates(.)[, "Y"]
    )
  
  # Verificación de variabilidad
  if (all(indice_parametro_i$valor == 0) || sd(indice_parametro_i$valor) < .Machine$double.eps) {
    # Crear un objeto que imite la estructura del resultado de autoKrige
    kriging_result <- list(
      krige_output = grilla_regular %>% dplyr::mutate(var1.pred = 0, var1.var = 0),
      exp_var = NULL,
      var_model = NULL
    )
    return(kriging_result)
  }
  
  # Interpolar parámetros usando kriging
  kriging_result <- automap::autoKrige(valor ~ 1, indice_parametro_i, grilla_regular)
  return(kriging_result)
}



EstimarParametrosConfiguracion <- function(input.value, 
                                           indice_parametro,
                                           ubicaciones_a_procesar,
                                           estaciones, 
                                           script,
                                           config) {
  
  # Obtener la configuracion para la cual se estimarán los parametros
  configuracion <- input.value
  
  # Convertir dataframe de estaciones en un objetp espacial 
  # para extraer las coordenadas de las estaciones
  estaciones.extraccion <- ubicaciones_a_procesar %>%
    sf::st_as_sf(., coords = c("lon_dec", "lat_dec"), crs = config$proj4string$latlon) %>%
    sf::st_transform(config$proj4string$planar)
  
  # Obtener parametros a estimar
  parametros <- indice_parametro %>%
    dplyr::filter(indice_configuracion_id %in% configuracion$id) %>%
    dplyr::pull(parametro) %>%
    unique()
  
  # Pentadas son parametros
  pentada_fin <- indice_parametro %>%
    dplyr::filter(indice_configuracion_id %in% configuracion$id) %>%
    dplyr::pull(pentada_fin) %>%
    unique()
  # Crear una grilla regular vacia para generar los mapas
  # a) Buscar zona de estudio regional y local
  area <- sf::read_sf(glue::glue("{config$dir$data}/shapefile/{config$files$area}")) %>%
    sf::st_transform(., config$proj4string$planar)  
  
  grilla_regular <- sf::st_bbox(area) %>%
    # Por defecto se usa una grilla regular de 10 km para la estimacion
    stars::st_as_stars(dx = 10000) %>%
    sf::st_crop(area) 
  
  # Iterar sobre las combinaciones de índices de configuración
  parametros_interpolados <- purrr::map(
    .x = configuracion$id,
    .f = function(indice) {
      
      script$info(glue::glue("Interpolando parametros para indice: {indice}\n"))
      
      # Iterar sobre las pentadas
      pentada_estimada <- purrr::map(
        .x = pentada_fin,
        .f = function(pentada) {
          
          script$info(glue::glue("Interpolando parametros para indice: {indice} y péntada: {pentada}.\n"))
          
          # Iterar sobre los parámetros
          parametro_estimados <- map(
            .x = parametros,
            .f = function(param) {
              
              script$info(glue::glue("Interpolando parametros para indice: {indice}, péntada: {pentada} y parámetro: {param}.\n"))
              
              # Estimar los parámetros
              resultado_kriging <- EstimarParametrosIndices(indice, pentada, param, estaciones, indice_parametro, grilla_regular, config)
              
              return(resultado_kriging)
            }
          ) %>% set_names(parametros) # Nombrar los elementos según los parámetros
          
          return(parametro_estimados)
        }
      ) %>% set_names(pentada_fin) # Nombrar los elementos según las pentadas
      
      return(pentada_estimada)
    }
  ) %>% set_names(unique(configuracion$internal_id)) # Nombrar los elementos según los índices de configuración
  
  # Iterar sobre los indices para guardar los parametros y rasters 
  # correspondientes a cada estacion
  purrr::map(
    .x = as.character(names(parametros_interpolados)),
    .f = function(indice) {
      
      indice.i <- parametros_interpolados[[indice]]
      
      purrr::map(
        .x = as.character(names(indice.i)),
        .f = function(pentada) {
          
          indice.pentada.i <- indice.i[[pentada]]
          
          parametros_pentada_ubicaciones <- purrr::map_dfr(
            .x = as.character(names(indice.pentada.i)),
            .f = function(param) {
              
              # Parametros con los rasters de cada 
              raster_parametros <- indice.pentada.i[[param]]$krige_output %>%
                terra::rast()
              names(raster_parametros) <- c("prediction", "var", "stdev")
              
              if (config$params$guardar.raster) {
                # Guardar resultados en un archivo
                results_filename <- glue::glue("{config$dir$data}/{config$files$indices_sequia$rasters}")
                results_filename <- stringr::str_replace(results_filename, "X", glue::glue("{indice}_{pentada}_{param}"))
                script$info(glue::glue("Guardando resultados en el archivo {results_filename}"))
                terra::writeRaster(raster_parametros, results_filename, overwrite = TRUE)
              }
              
              parametros_estimados_ubicaciones <- raster_parametros %>%
                terra::extract(estaciones.extraccion) %>%
                dplyr::select(-ID) %>%
                dplyr::mutate(station_id = ubicaciones_a_procesar$station_id,
                              parametro = param) %>%
                tidyr::gather(variable, valor, -station_id, -parametro)
              
            }
          )
          
          # Iterar por todas las estaciones para guardar los parametros
          for (station in unique(ubicaciones_a_procesar$station_id)) {
            
            # Guardar resultados en un archivo
            results_filename <- glue::glue("{config$dir$data}/{config$files$indices_sequia$parametros_intermedios}_{station}_{indice}_{pentada}.feather")
            script$info(glue::glue("Guardando resultados en el archivo {results_filename}"))
            
            # Extraer el numero de configuracion correspondiente al indice
            id_configuracion <- configuracion %>%
              dplyr::filter(internal_id == !!indice) %>%
              dplyr::pull(id)
            
            parametros_pentada_ubicaciones %>%
              dplyr::filter(station_id == !!station,
                            variable == "prediction") %>%
              dplyr::mutate(metodo_imputacion = 0,
                            pentada.fin = pentada,
                            indice_configuracion_id = id_configuracion) %>%
              # Corregir probabilidad de 0
              dplyr::mutate(valor = if_else(parametro == "prob.0" & valor < 0, 
                                            0 , valor)) %>%
              dplyr::mutate(valor = if_else(parametro == "prob.media.0" & valor < 0, 
                                            0 , valor)) %>%
              feather::write_feather(., results_filename)
            
            # Si se selecciona la simulacion de parametros para estimar la variabilidad
            # de los parametros se realizarán n simulaciones de una distribucion normal
            # considerando los valores interpolados y sus errores respectivos de estimacion
            if (config$params$simulacion.parametros$simular) {
              
              script$info(glue::glue("Simulando parametros a partir de la estaciones para la estacion: {station} y pentada: {pentada}"))
              
              # Fijar una semilla para asegurar reproducibilidad
              set.seed(1234)
              
              parametros_simulados <- parametros_pentada_ubicaciones %>%
                dplyr::filter(station_id == !!station,
                              variable != "var") %>%
                tidyr::pivot_wider(names_from = variable, values_from = valor) %>%
                dplyr::rowwise() %>%
                dplyr::mutate(simulacion = list(rnorm(config$params$simulacion.parametros$n_muestras,
                                                      mean = prediction, sd = stdev))) %>%
                tidyr::unnest(cols = c(simulacion)) %>%
                dplyr::group_by(parametro) %>%
                dplyr::mutate(id_realizacion = 1:config$params$simulacion.parametros$n_muestras) %>%
                dplyr::select(station_id, parametro, id_realizacion, everything()) %>%
                dplyr::ungroup()
              
              # Simular N muestras de paramentros para la estacion elegida
              for (realizacion in unique(parametros_simulados$id_realizacion)) {
                
                # Guardar resultados en un archivo
                results_filename <- glue::glue("{config$dir$data}/{config$files$indices_sequia$parametros_intermedios}_{station}_{indice}_{pentada}_{realizacion}.feather")
                script$info(glue::glue("Guardando resultados en el archivo {results_filename}"))
                
                # Extraer el numero de configuracion correspondiente al indice
                id_configuracion <- configuracion %>%
                  dplyr::filter(internal_id == !!indice) %>%
                  dplyr::pull(id)
                
                parametros_simulados %>%
                  dplyr::filter(station_id == !!station,
                                id_realizacion == !!realizacion) %>%
                  # Corregir probabilidad de 0
                  dplyr::mutate(simulacion = if_else(parametro == "prob.0" & simulacion < 0, 
                                                     0 , simulacion)) %>%
                  dplyr::mutate(simulacion = if_else(parametro == "prob.media.0" & simulacion < 0, 
                                                     0 , simulacion)) %>%
                  dplyr::select(station_id, parametro, valor = simulacion) %>%
                  dplyr::mutate(metodo_imputacion = 0,
                                pentada.fin = pentada,
                                indice_configuracion_id = id_configuracion) %>%
                  feather::write_feather(., results_filename)
              }
              
            } 
            
          }
        }
      )
    }
  )
}
# ------------------------------------------------------------------------------