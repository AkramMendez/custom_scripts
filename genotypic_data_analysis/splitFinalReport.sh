#!/bin/bash
#
# Developer:  Akram Méndez
# Date:       17 Jul 2018
# Script:     Extractor de datos para cada muestra a partir de un FinalReport.
#             Código 46, 2018.
#

final_reports=$1
output_path=$2
type=$3

[ $# -eq 0 ] && { echo "
                Extractor de variantes para cada muestra.
                Script: 
                                      Este script se utiliza para generar un archivo para cada muestra contenida en un lote de FinalReports.
                                      Cada archivo generado contiene todas las variantes del chip (custom v1).

                Uso: 
                                      $0 <raw_FinalReport o raw_FinalReports_dir > <output_path (optional)> <file_type: FZERO or empty for default (optional)>
                
                Definición: 
                FinalReport(s)          Ruta al archivo o directorio que contiene los final-reports crudos provenientes de GenomeStudio.
                Project_path (opcional) Ruta del directorio donde se guardarán los archivos generados. 
                                        La ruta por default al omitir este parámetro es:
                                        /media/datashare/data/Produccion/reports-all-variants

                Ejemplo: 
                                      ./getSampleData.sh final-reports-raw/final-report.txt output_path FZERO
 "; exit 1;}

if [ -z "$2" ]
  then
  echo "Setting output path to $output_path"
  output_path="/media/datashare/data/Produccion/reports-all-variants"
fi;

if [ -z "$3" ]
  then
  echo "Setting file type to default format"
  type="FZERO"
fi;

splitFinalReport(){
	local report=$1
  local output=$2
	awk -v outpath="$output" 'BEGIN{FS="\t";OFS="\t"}{if(NR<=10){next} else{ print "-",NR,$2; print $0 >> outpath"/"$2"_all_variants.txt"}}' $report
}


filename=$(realpath ${final_reports})
project_name=$( echo "${filename##*/}" | sed -E 's/(_FinalReport.*.txt|.txt)//g')

logfile="/media/datashare/data/Produccion/logfiles/logfile_splitted_FinalReports.txt"

path=$( echo "${output_path}/${project_name}")

if [ type=="FZERO" ]
  then
  header=$(echo -e "SNP Name\tSample ID\tAllele1 - Top\tAllele2 - Top\tAllele1 - Forward\tAllele2 - Forward\tGC Score\tB Allele Freq\tLog R Ratio")
else
  header=$(echo -e "SNP Name\tSample ID\tAllele1 - Top\tAllele2 - Top\tAllele1 - Forward\tAllele2 - Forward\tGC Score")
fi;

# Check if the output directory exists
if [ ! -d "${path}" ]   
  then
    echo "Making dir: ${path}"
    mkdir -p "${path}"
  fi;

# Extract data for a single FinalReport or for multiple FinalReports in a directory
if [ -f ${final_reports} ]
  then
echo "Processing FinalReport: ${final_reports}"
echo "Converting file to unix format (dos2unix):"
dos2unix $final_reports
echo "Splitting FinalReport"
splitFinalReport $final_reports $path
fi; 

if [ -d ${final_reports} ]
  then
echo "Processing FinalReports in: ${final_reports}"
    for i in $(ls ${final_reports}); do
      echo "Converting file to unix format (dos2unix):"
      dos2unix $i
      echo "Splitting FinalReport"
      splitFinalReport $i $path
    done
fi; 

echo "Appending header to FinalReports"
#Add header to the FinalReport's first line
  find "${path}" -type f -exec sed -i "1s/^/${header}\n/" {} \;  
echo "FinalReport saved to: $path"
echo "Saving file name to logfile"
echo -e "${project_name}\t$(date +'%Y%m%d')" >> ${logfile}  
echo "Done."
