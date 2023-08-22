#!/usr/bin/env bash
# Author:       Akram Méndez
# Email:        amendez@codigo46.com.mx
# Date:         2020/03/17
# Copyright: Codigo46 S.A. de C.V.
#


[ $# -eq 0 ] && { echo "
                Compilador de codigos de Extracción
                Script:
                                      Este programa se utiliza para recopilar nuevos lotes de extracción añadiéndolos a la bitácora de extracciones general 
                                      que contiene todas las muestras procesadas hasta el momento.
                                      
                Uso:
                                      $0 <archivo de Extracción>

                Definición:
                archivo de Extracción     Archivo del lote de Extracción más reciente, generado en el laboratorio, previamente sincronizado a la instancia de producción (/media/datashare/data/Produccion/extractions)
                
                Ejemplo:
                                      compile_extractions.sh Extraccion_20200220_CVS_LAC.csv
 "; exit 1;}


extraction_file=$1
extractions_compilation=$(realpath /media/datashare/data/Produccion/extractions/*compilado-extracciones-al-*)
compilation_backup='/media/datashare/data/Produccion/extractions/testing/all_extractions_updated.txt'

echo "Converting extraction file (dos2unix): ${extraction_file}"
dos2unix ${extraction_file}

#Obtener nueva fecha para actualizar los compilados de extracción (se toma la fecha del nombre del útimo archivo de extracción):
date=$(ls $1 | grep -oP "(?<=\_)(\d+)(?=\_)")
new_file=$(ls ${extractions_compilation} | grep -oP "(.*\-)(?=\d+\.txt)" )
output=$(echo "${new_file}${date}.txt")

# Check if the output directory exists
if [ ! -f "${output}" ] 
then
  echo "Actualizando compilado de Extracción:"
  echo ${extraction_file}
  awk 'BEGIN{FS=",|\t|;";OFS="\t"}NR>1{print $2,$3}' ${extraction_file} >> ${extractions_compilation}
  echo "Actualizando backup del compilado de Extracción:"
  echo ${compilation_backup}
  awk 'BEGIN{FS=",|\t|;";OFS="\t"}NR>1{print $2,$3}' ${extraction_file} >> ${compilation_backup}
  echo "Actualizando nombre del archivo a la fecha de la última extracción."
  mv ${extractions_compilation} ${output}
else
  echo "El compilado de extracción ya está actualizado a la última extracción."
fi


