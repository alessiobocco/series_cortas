# -----------------------------------------------------------------------------#
# --- PASO 1. Cargar paquetes necesarios ----
rm(list = ls()); gc()
Sys.setenv(TZ = "UTC")
list.of.packages <- c("ADGofTest", "caret", "dplyr", "fitdistrplus", "lmomco", "logspline",
                      "lubridate", "magrittr", "mgcv", "purrr", "SCI", "sirad", "SPEI", 
                      "stats", "stringr", "utils", "WRS2", "yaml", "yardstick", "feather",
                      "R6", "futile.logger", "mgcv", "doSNOW", "foreach", "snow", "parallel",
                      "RPostgres", "data.table", "purrr", "automap", "sf")
for (pack in list.of.packages) {
  if (!require(pack, character.only = TRUE)) {
    stop(paste0("Paquete no encontrado: ", pack))
  }
}
rm(pack); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- PASO 2. Leer archivo de configuracion ----
# -----------------------------------------------------------------------------#

normalize_dirnames <- function(dirnames) {
  if (is.atomic(dirnames)) 
    dirnames <- base::sub('/$', '', dirnames)
  if (!is.atomic(dirnames))
    for (nm in names(dirnames)) 
      dirnames[[nm]] <- normalize_dirnames(dirnames[[nm]])
  return (dirnames)
}

# a) YAML de configuracion del cálculo de índices de sequia
args <- base::commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  archivo.config <- args[1]
} else {
  # No vino el archivo de configuracion por linea de comandos. Utilizo un archivo default
  archivo.config <- paste0(getwd(), "/configuracion_calculador_indices.yml")
}
if (! file.exists(archivo.config)) {
  stop(paste0("El archivo de configuración ", archivo.config, " no existe\n"))
} else {
  cat(paste0("Leyendo archivo de configuración ", archivo.config, "...\n"))
  config <- yaml::yaml.load_file(archivo.config)
  config$dir <- normalize_dirnames(config$dir)
  
}

# b) YAML de parametros del cálculo de índices de sequia
if (length(args) > 1) {
  archivo.params <- args[2]
} else {
  # No vino el archivo de configuracion por linea de comandos. Utilizo un archivo default
  archivo.params <- paste0(getwd(), "/parametros_calculador_indices.yml")
}
if (! file.exists(archivo.params)) {
  stop(paste0("El archivo de parámetros ", archivo.params, " no existe\n"))
} else {
  cat(paste0("Leyendo archivo de parámetros ", archivo.params, "...\n"))
  config$params <- yaml::yaml.load_file(archivo.params)
}

# c) YAML de configuración de la API
if (length(args) > 1) {
  archivo.config.api <- args[3]
} else {
  # No vino el archivo de configuracion por linea de comandos. Utilizo un archivo default
  archivo.config.api <- paste0(getwd(), "/configuracion_api.yml")
}
if (! file.exists(archivo.config.api)) {
  stop(paste0("El archivo de configuración ", archivo.config.api, " no existe\n"))
} else {
  cat(paste0("Leyendo archivo de configuración ", archivo.config.api, "...\n"))
  config$api <- yaml::yaml.load_file(archivo.config.api)$api
}

replace_run_identifier <- function(filenames, identifier) {
  if (is.atomic(filenames) && grepl("<\\*idc>", filenames)) 
    filenames <- base::sub('<\\*idc>', identifier, filenames)
  if (!is.atomic(filenames))
    for (nm in names(filenames)) 
      filenames[[nm]] <- replace_run_identifier(filenames[[nm]], identifier)
  return (filenames)
}

# d) YAML de configuración del intercambio de archivos del proceso de generación de índices
if (length(args) > 1) {
  archivo.nombres <- args[3]
} else {
  # No vino el archivo de configuracion por linea de comandos. Utilizo un archivo default
  archivo.nombres <- paste0(config$dir$data, "/configuracion_archivos_utilizados.yml")
}
if (! file.exists(archivo.nombres)) {
  stop(paste0("El archivo de configuración ", archivo.nombres, " no existe\n"))
} else {
  cat(paste0("Leyendo archivo de configuración ", archivo.nombres, "...\n"))
  config$files <- yaml::yaml.load_file(archivo.nombres)
  config$files <- replace_run_identifier(config$files, config$files$identificador_corrida)
}

rm(archivo.config, archivo.params, archivo.nombres, args); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- PASO 3. Cargar librerias propias e iniciar script ----
# -----------------------------------------------------------------------------#

# a) Cargar librerias
source(glue::glue("{config$dir$lib}/FechaUtils.R"), echo = FALSE)
source(glue::glue("{config$dir$lib}/Script.R"), echo = FALSE)
source(glue::glue("{config$dir$lib}/Task.R"), echo = FALSE)
source(glue::glue("{config$dir$lib}/Helpers.R"), echo = FALSE)
source(glue::glue("{config$dir$lib}/crc-api.R"), echo = FALSE)
source(glue::glue("{config$dir$lib}/Servicios.R"), echo = FALSE)

# b. Carga de codigo para ajuste de distribuciones, calculo de indices y ejecucion distribuida
source(glue::glue("{config$dir$base}/lib/funciones_ajuste.R"), echo = FALSE)
source(glue::glue("{config$dir$base}/lib/funciones_calculo.R"), echo = FALSE)
source(glue::glue("{config$dir$base}/lib/funciones_test_ajuste.R"), echo = FALSE)
source(glue::glue("{config$dir$base}/lib/funciones_worker.R"), echo = FALSE)

# c) Chequear que no este corriendo el script de estadisticas.
#    Si esta corriendo, la ejecucion debe cancelarse.
script.estadisticas <- Script$new(run.dir = config$dir$estadisticas$run,
                                   name = "EstadisticaMovil")
script.estadisticas$assertNotRunning()
rm(script.estadisticas)


# b.1) Definir nombre del script
script_name <- "IndicesGenerador"
script_logfile <- glue::glue("{config$dir$run}/{script_name}.log")

# b.2) borrar archivo .log de corridas anteriores
if (file.exists(script_logfile))
  file.remove(script_logfile)

# b.3) Iniciar script
# Si el directorio run para almacenar los log no existe, crearlo
if (!fs::dir_exists(config$dir$run)) {
  fs::dir_create(config$dir$run)
}

# Si el directorio para almacenar los resultados de los 
# tests de bondad de ajuste, crearlo
if (!fs::dir_exists(glue::glue("{config$dir$data}/control/resultados_tests"))) {
  fs::dir_create(glue::glue("{config$dir$data}/control/resultados_tests"))
}
# Si el directorio para almacenar los parametros de las 
# funciones de distribución de probabilidad, crearlo
if (!fs::dir_exists(glue::glue("{config$dir$data}/control/parametros"))) {
  fs::dir_create(glue::glue("{config$dir$data}/control/parametros"))
}
script <- Script$new(run.dir = config$dir$run, name = script_name, create.appender = T)
script$start()

# Crear directorios para guardar resultados intermedios
if (!fs::dir_exists(glue::glue("{config$dir$data}/partial"))) {
  fs::dir_create(glue::glue("{config$dir$data}/partial"))
}

# Crear directorios para guardar resultados intermedios: parametros
if (!fs::dir_exists(glue::glue("{config$dir$data}/control/parametros"))) {
  fs::dir_create(glue::glue("{config$dir$data}/control/parametros"))
}

# Crear directorios para guardar resultados intermedios: resultados_tests
if (!fs::dir_exists(glue::glue("{config$dir$data}/control/resultados_tests"))) {
  fs::dir_create(glue::glue("{config$dir$data}/control/resultados_tests"))
}

# Crear directorios para guardar resultados finales
if (!fs::dir_exists(glue::glue("{config$dir$data}/output"))) {
  fs::dir_create(glue::glue("{config$dir$data}/output"))
}

# Crear directorios para guardar resultados finales
if (!fs::dir_exists(glue::glue("{config$dir$data}/output/rasters"))) {
  fs::dir_create(glue::glue("{config$dir$data}/output/rasters"))
}
# e) Cargar funciones necesarias
source(glue::glue("../../lib/crc-api.R"), echo = FALSE)

# f) Variable para almacenar los posibles errores
errors <- c()

# Configuracion SSL
httr::set_config( httr::config(ssl_verifypeer = FALSE))
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- PASO 4. Inicializar variables ----
# -----------------------------------------------------------------------------#

# a) Leer datos necesario de la BBDD
httr::set_config( httr::config(ssl_verifypeer = FALSE) )

# Obtenemos los datos de las estaciones.
estaciones <- purrr::map_dfr(
  .x = config$params$paises,
  .f = function(pais) {
    est <- ConsumirServicioJSON(url = paste0(config$api$url, glue::glue("/estaciones/{pais}")),
                                usuario = config$api$user, clave = config$api$pass) 
    est <- est %>%
      dplyr::mutate(pais_id = pais) 
    return (est)
  }
) 

# Filtrar solo las estaciones convencionales
estaciones %<>% 
  dplyr::filter(tipo == "C") %>%
  dplyr::rename(lat_dec = latitud, lon_dec = longitud)

# Configuraciones de indices
#indice_configuracion <- dplyr::tbl(con, "indice_configuracion") %>%
#  dplyr::collect() %>%
#  dplyr::rename(internal_id = id) %>%
#  dplyr::mutate(indice = as.character(indice),
#                distribucion = as.character(distribucion),
#                metodo_ajuste = as.character(metodo_ajuste))
# readr::write_csv(indice_configuracion, "indice_configuracion.csv")

# Fitrar configuraciones deseadas
indice_configuracion <- readr::read_csv("./indice_configuracion.csv") %>%
  dplyr::filter(#escala %in% config$params$escalas,
                escala %in% 1,
                indice %in% "SPI",
                metodo_ajuste == "ML-SinRemuestreo")

# Crear objeto con los configuraciones para cada estacion e indice
combinaciones.parametros <- tidyr::expand_grid(omm_id = estaciones$omm_id,
                   configuracion_id = indice_configuracion$internal_id)

# Consultar parámetros para las estaciones seleccionadas
indice_parametro <- purrr::map2_dfr(
  .x = combinaciones.parametros$omm_id,
  .y = combinaciones.parametros$configuracion_id,
  .f = function(estacion_id, configuracion_id) {
    
    parametros <- ConsumirServicioJSON(url = paste0(config$api$url, glue::glue("/indices_sequia_parametros_ajuste/{configuracion_id}/{estacion_id}")),
                                usuario = config$api$user, clave = config$api$pass) 
    
  }
)

# c) Buscar las estadisticas moviles 
script$info(glue::glue("Buscando estadísticas móviles para calcular indices de sequia, ",
                       "archivo: {config$files$estadisticas_moviles$resultados}"))
archivo <- glue::glue("{config$dir$data}/{config$files$estadisticas_moviles$resultados}")
estadisticas.moviles <- data.table::fread(archivo, nThread = config$files$avbl_cores) %>%
  dplyr::mutate(fecha_desde = as.Date(fecha_desde), fecha_hasta = as.Date(fecha_hasta)) %>%
  tibble::as_tibble()
base::remove(archivo); base::invisible(base::gc())
script$info("Seleccionando estadísticas móviles asociadas a las escalas especifícadas")
estadisticas.moviles <- estadisticas.moviles %>%
  dplyr::filter(ancho_ventana_pentadas %in% (config$params$escalas * 6))

# d) Verificar que hayan estadísticas para todas las escalas especificadas en el yaml
if (!all((config$params$escalas * 6) %in% estadisticas.moviles$ancho_ventana_pentadas))
  stop("En parametros_calcular_indices.yml se han especificado escalas para las ",
       "cuales no han sido calculadas estadísticas móviles!!")

# e) Obtener configuraciones para el cálculo de los indices de sequía
script$info(glue::glue("Buscando configuraciones para los índices a ser calculados, ",
                       "archivo: {config$files$indices_sequia$configuraciones}"))
archivo <- glue::glue("{config$dir$data}/{config$files$indices_sequia$configuraciones}")
configuraciones.indices <- data.table::fread(archivo, nThread = config$files$avbl_cores) %>% 
  dplyr::mutate(referencia_comienzo = as.Date(referencia_comienzo), referencia_fin = as.Date(referencia_fin)) %>%
  tibble::as_tibble()
base::remove(archivo); base::invisible(base::gc())
script$info("Seleccionando configuraciones de índice asociadas a las escalas especifícadas")
configuraciones.indices <- configuraciones.indices %>%
  dplyr::filter(escala %in% config$params$escalas)

# f) Verificar que hayan configuraciones para todas las escalas especificadas en el yaml
configuraciones.indices <- configuraciones.indices %>%
  dplyr::group_by(indice, distribucion, metodo_ajuste) %>%
  dplyr::group_walk(.f = function(g, k) {
      if (!all(config$params$escalas %in% g$escala)) {
        stop_msg <- glue::glue("Algunas de las escalas definidas en el archivo parametros_calculador.indices.yml ",
                               "no están presentes para la configuración con indice:{k$indice}, distribucion:",
                               "{k$distribucion} y metodo_ajuste:{k$metodo_ajuste}!!")
        stop(stop_msg)
      }
    })

# Agregar indice interno para los indices de sequia
configuraciones.indices %<>%
  dplyr::left_join(indice_configuracion) 
rm(indice_configuracion)

# g) Buscar ubicaciones a las cuales se aplicara el calculo de indices de sequia
# g.1) Obtener datos producidos por el generador y filtrarlos
script$info("Leyendo csv con datos de entrada")
csv_filename <- glue::glue("{config$dir$data}/{config$files$clima_generado}")
datos_climaticos <- data.table::fread(csv_filename, nThread = config$files$avbl_cores) 
script$info("Lectura del csv finalizada")
# g.x) Reducción de trabajo (solo para pruebas)
# datos_climaticos_generados <- datos_climaticos_generados %>%
#   dplyr::filter( realization %in% c(1, 2, 3), dplyr::between(date, as.Date('1991-01-01'), as.Date('2000-12-31')) )
# g.2) Generar tibble con ubicaciones sobre las cuales iterar
script$info("Obtener ubicaciones sobre las cuales iterar")
ubicaciones_a_procesar <- datos_climaticos %>%
  dplyr::select(dplyr::ends_with("_id"), 
                lon_dec = longitude, 
                lat_dec = latitude) %>%
  tibble::as_tibble() %>% dplyr::distinct()
script$info("Obtención de ubicaciones a iterar finalizada")
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- PASO 4. Estimar parametros de manera distriuída ----
# -----------------------------------------------------------------------------#
function_name <- "EstimarParametrosConfiguracion"
task_logfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.log")
task_outfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.out")

# Borrar archivos .log y .out de corridas anteriores
if (file.exists(task_logfile))
  file.remove(task_logfile)
if (file.exists(task_outfile))
  file.remove(task_outfile)

# Crear tarea distribuida y ejecutarla
task.estimar.parametros <- Task$new(parent.script = script,
                                func.name = function_name,
                                packages = list.of.packages)

# Ejecutar tarea distribuida
script$info("Estimando parámetros de indices de sequia")
resultados.indices.sequia <- task.estimar.parametros$run(number.of.processes = config$max.procesos, 
                                                     config = config, 
                                                     script = script,
                                                     estaciones = estaciones, 
                                                     ubicaciones_a_procesar = ubicaciones_a_procesar,
                                                     #input.values = configuraciones.indices,
                                                     input.value = configuraciones.indices[1,],
                                                     indice_parametro = indice_parametro)
# Agregar log de la tarea al log del script
file.append(script_logfile, task_logfile)

# Si hay errores, terminar ejecucion
task.estimar.parametros.errors <- task.estimar.parametros$getErrors()
if (length(task.estimar.parametros.errors) > 0) {
  for (error.obj in task.estimar.parametros.errors) {
    id_column <- IdentificarIdColumn(ubicaciones_a_procesar %>% dplyr::top_n(1))
    script$warn(glue::glue("({id_column}={error.obj$input.value[[id_column]]}): {error.obj$error}"))
  }
  script$error("Finalizando script de forma ANORMAL")
}

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- PASO 5. Calcular indices de sequia de manera distribuida ----
# -----------------------------------------------------------------------------#

# Definir nombre de la función a ser distribuida y nombre de archivos .log y .out de corridas anteriores
function_name <- "CalcularIndicesSequiaUbicacion"
task_logfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.log")
task_outfile <- glue::glue("{config$dir$run}/{script_name}-{function_name}.out")

# Borrar archivos .log y .out de corridas anteriores
if (file.exists(task_logfile))
  file.remove(task_logfile)
if (file.exists(task_outfile))
  file.remove(task_outfile)

# Crear tarea distribuida y ejecutarla
task.indices.sequia <- Task$new(parent.script = script,
                                func.name = function_name,
                                packages = list.of.packages)

# Ejecutar tarea distribuida
script$info("Calculando indices de sequia")
resultados.indices.sequia <- task.indices.sequia$run(number.of.processes = config$max.procesos, 
                                                     config = config, 
                                                     input.values = ubicaciones_a_procesar, 
                                                     #input.value = ubicaciones_a_procesar[1,],
                                                     configuraciones.indices, estadisticas.moviles)

# Agregar log de la tarea al log del script
file.append(script_logfile, task_logfile)

# Si hay errores, terminar ejecucion
task.indices.sequia.errors <- task.indices.sequia$getErrors()
if (length(task.indices.sequia.errors) > 0) {
  for (error.obj in task.indices.sequia.errors) {
    id_column <- IdentificarIdColumn(ubicaciones_a_procesar %>% dplyr::top_n(1))
    script$warn(glue::glue("({id_column}={error.obj$input.value[[id_column]]}): {error.obj$error}"))
  }
  script$error("Finalizando script de forma ANORMAL")
}

# Transformar resultados a un objeto de tipo tibble
resultados.indices.sequia.tibble <- resultados.indices.sequia %>% purrr::map_dfr(~.x)

# Guardar resultados en un archivo fácil de compartir
results_filename <- glue::glue("{config$dir$data}/{config$files$indices_sequia$resultados}")
script$info(glue::glue("Guardando resultados en el archivo {results_filename}"))
data.table::fwrite(resultados.indices.sequia.tibble, file = results_filename, nThread = config$files$avbl_cores)

# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- PASO 6. Finalizar script cerrando conexion a base de datos ----
# -----------------------------------------------------------------------------#

# a) Finalizar script
script$stop()

# b) Crear archivo .info
info_filename <- glue::glue("{config$dir$data}/{config$files$indices_sequia$info_corrida}")
if (file.exists(info_filename))
  file.remove(info_filename)
file.copy(from = script_logfile, to = info_filename)

# ------------------------------------------------------------------------------
