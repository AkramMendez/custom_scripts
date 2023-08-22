#!/bin/bash
#
# Developer:  Akram Méndez
# Date:       28 Ene 2019
# Script:     Generador de archivos "chipMuestra" a partir de Samples_Tables para lotes de muestras Infinium.
#             Código 46, 2019.


file=$1


[ $# -eq 0 ] && { echo "
                Generador de archivos 'chipMuestra' a partir de archivos Samples_Tables.
                Script: Este script genera un archivo 'chipMuestra' a partir de los Sample_IDs, Sentrix_ID y Sentrix_positions
                		obtenidos a partir de un archivo de Samples_Tables correspondiente a un proyecto Infinium.
                		El archivo chipMuestra se usa posteriormente en la base de datos para actualizar la información de los chips y las muestras relacionadas.
                Uso: 
                                      $0 <Samples_Table>
                
                Definición: 
                <Samples_Table>		Tabla 'Samples_Table' con información de muestras, chip, CR scores, etc. exportada a partir de un proyecto Infinium

                Ejemplo: 
                                      ./makeChipMuestra Samples_Table_Infinium_20190122_BTH_ERM.txt
 "; exit 1;}


chipMuestra(){
	awk 'BEGIN{FS="\t";OFS=","}{if(NR==1){gsub("\r$", "");print "SAMPLE_ID","Array Info.Sentrix ID","Array Info.Sentrix Position"}else{print $6,$24,$25}}' ${1}
}


if [ -f ${file} ]
  then
  chipMuestra ${file}
fi; 

if [ -d ${file} ]
  then
  for i in $(ls ${file}); do
  chipMuestra ${i} ;done    
fi;

