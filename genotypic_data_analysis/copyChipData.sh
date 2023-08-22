#!/usr/bin/env bash
# Author:       Akram Méndez
# Email:        amendez@codigo46.com.mx
# Date:         2020/03/17
# Copyright: Codigo46 S.A. de C.V.
#


[ $# -eq 0 ] && { echo "
                Copiador de datos del LIMS para muestras de Usuarios por Volumen.
                Script:
                                      Este programa se utiliza para buscar y copiar los archivos de entrega para Usuarios por Volumen a partir de una lista de códigos 
                                      'sentrix_id' y 'sentrix_position' correspondientes a las placas de procesadas para dicho lote.
                                      Para cada muestra con un sentrix_id y sentrix_position a partir del archivo Samples_Map, 
                                      el programa crea las siguientes carpetas en la ruta de salida y recopila la infomación correspondiente:

                                      1. '/decode/': 
                                      2. '/image_data/':
                                      3. '/intensity_data':
                Uso:
                                      $0 <Archivo samples_Map (TSV)> <Illumina LIMS mounting path> <Output path>

                Definición:
                samples_Map          Archivo TSV con la siguiente información: [Columna1 = 'Sample ID' ;Columna2 = 'sentrix_id'; Columna3 = 'sentrix_position']
                Illumina LIMS path   Ruta a la partición donde está montado el LIMS (A partir de esta ruta buscará las carpetas 'decode','image_data' e 'intensity_data')
                output_path          Ruta del proyecto donde se almacenarán los resultados de cada lote (Usualmente se utiliza la carpeta 'bioinfobatch' en el LIMS)

                Ejemplo:
                                      ./copyChipData.sh Samples_Map.txt /media/sambaLims/ /media/sambaLims/bioinfobatch/muestras_INMEGEN/
 "; exit 1;}

samples_map=$1
lims_path=$(realpath "$2")
output_path=$(realpath "$3")

sentrix_ids=$(awk 'BEGIN{FS=",|\t|;";OFS="\t"}NR>1{print $2}' ${samples_map} | sort | uniq)

decode="decode"
image_data="image_data"
intensity_data="intensity_data"

# Check if the output directory exists
directoryExists(){
  if [ ! -d "$1" ]
  then
    echo "Making dir: $1"
    mkdir -p "$1"
  fi;
}

copyFolder(){
  find "$1/$2" -iname "${3}_${4}*" -exec cp {} $5 \;
}


echo "Checking output directories ..."
directoryExists ${decode}
directoryExists ${intensity_data}
directoryExists ${image_data}

for sentrix_id in ${sentrix_ids[@]}; do
  echo "Copying sentrix_id: ${sentrix_id}"
 directoryExists "${decode}/${sentrix_id}"
  directoryExists "${image_data}/${sentrix_id}"
  directoryExists "${intensity_data}/${sentrix_id}"
  for sentrix_position in $(cat ${samples_map}  | grep ${sentrix_id} | cut -d, -f3 | sort); do
    echo "Copying sentrix_position: ${sentrix_position}"
    copyFolder ${lims_path} "decode" ${sentrix_id} ${sentrix_position} "${output_path}/${decode}/${sentrix_id}"
    copyFolder ${lims_path} "image_data" ${sentrix_id} ${sentrix_position} "${output_path}/${image_data}/${sentrix_id}"
    copyFolder ${lims_path} "intensity_data" ${sentrix_id} ${sentrix_position} "${output_path}/${intensity_data}/${sentrix_id}";done;done

    echo "Done."
