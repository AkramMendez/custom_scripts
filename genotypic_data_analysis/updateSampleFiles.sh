#!/bin/bash
#
# Developer:  Akram Méndez
# Date:       14 Ene 2019
# Script:      Generador de archivos 'updateSamples' a partir de archivos de Extracción.
#             Código 46, 2019.
#

extraction_file=$1


[ $# -eq 0 ] && { echo "
                Generador de archivos 'Samples' a partir de archivos de Extracción.
                Script: Este script genera un archivo 'Samples' a partir de los Sample_IDs y barCodes 
                		obtenidos a partir de un archivo de extracción concatenando dichos códigos en
                		el esqueleto de un query para la posterior actualización de las muestras en la 
                		base de datos. 
                		                                     
                Uso: 
                                      $0 <extractionFile>
                
                Definición: 
                <extraction_file>	Archivos de extracción conteniendo los códigos Sample_ID y barcodes para un lote de muestras (Localizados en la carpeta predeterminada ftp://codigo46.me/Laboratorio/Extraccion).

                Ejemplo: 
                                      ./updateSamplesFile.sh Extraction_20190101_AMR_RGR.csv
 "; exit 1;}


updateSamplesFile(){	
  awk 'BEGIN{FS="[,|;]";OFS=""}NR>1{gsub("\r$", "");print "UPDATE Samples SET Sample_ID = ","\x27",$2,"\x27"," WHERE BarCode = ","\x27",$3,"\x27"," AND Sample_ID = \x27\x27\n","UPDATE EstadoSamples SET FK01_Estado = 6 WHERE BarCode = ","\x27",$3,"\x27"'} ${1}
}



if [ -f ${extraction_file} ]
  then
 updateSamplesFile ${extraction_file}
fi; 

if [ -d ${extraction_file} ]
  then
  for i in $(ls ${extraction_file}); do
    updateSamplesFile ${i} ;done    
fi;

