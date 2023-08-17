
#!/bin/bash

#
# It makes KS analysis of the indicated pairs of categories, after MonteCarlo simulation of 
# t-Student distributions of all per-category replicas (technical and biological) to estimate the
# most likely DeltaDeltaCT values. At the end, 
#
# the script also includes a list of KS-tests across all possible pairs
#
# This script can be executed after make_qpcr_analysis1.sh
#
# USAGE
#
# bash make_qpcr_analysis2.sh --wd <string> -- analysis.folder <string> ...
#
# --wd <string> : This is the top working directory (absolute path) were global folders for data and results can be found. 
# --analysis.folder <string> : name of the analysis folder where results will be stored. Avoid spaces for compatibility!
# --data.folder <string> : name of the data folder. Compatible files are in CSV format, with tab as separator and quoted text. Decimal character is a dot. Avoid spacies for compatibility!
# --exp.folder <string>: name of the experimental folder where experimental file and reference files are located. Avoid spaces ...
# --exp.file <string>: name of the experimental file from QS3. Avoid spaces ...
# --ref.file <string>: name of the reference gene file from QS3. Avoid spaces ...
# --exp.levels <string> : name of the categories that will be normalised to the control in the same position (check --ctr.levels)
# --ctr.levels <string> : control names that will be compared to the specified exp.levels
# --montecarlo.n <integer> : number of points in the montecarlo simulation
# --qranks <integer,integer> : probability levels used to calculate the quantiles.
# --wgraph <float> : width of the graphical output
# --wh <float,float> : width and height in inches of the PDF output
# --lve <string> : composed levels used to make comparisons with the level at the same position in --cve
# --cve <string> : composed levels used to make comparisons with the level at the same position in --lve
# --dig <integer> : number of significant digits in the latex table
#
# @AUTHOR: Maximo Sanchez-Aragon
# @EMAIL: maxsanara@gmail.com
# 


declare -r RED_FORMAT="$(printf '\033[31;1m')"
declare -r BLUE_FORMAT="$(printf '\033[34;1m')"
declare -r FORMAT_RESET="$(printf '\033[0m')"
declare -r CONFIGFILE_NAME="config_ref"
declare -r SCRIPTPATH=$(dirname -- "$( readlink -f -- "$0"; )");
declare -r CONFIGFILE="${SCRIPTPATH}/${CONFIGFILE_NAME}"


# FUNCTIONS ----------------------
function read_property {
    PROPERTIES_FILE="${1}"
    while read -r line
    do
     PROP_NAME=$(echo $line | cut -f1 -d=)
     PROP_VALUE=$(echo $line | cut -f2 -d=)

     if [ $PROP_NAME == "${2}" ]; then
        echo $PROP_VALUE
     fi
   
    done < "$PROPERTIES_FILE"
}

echo "${BLUE_FORMAT}[INFO]${FORMAT_RESET} Reading properties file at ${CONFIGFILE} ..."

# ARGUMENTS assignation  ...
declare WD="${HOME}" #$(read_property "${CONFIGFILE}" "WD") 	# The working directory containing the analysis and data folders
declare ANFLD=$(read_property "${CONFIGFILE}" "ANFLD")	# The analysis folder name contained in ${WD}
declare DATAFLD=$(read_property "${CONFIGFILE}" "DATAFLD")  # The data folder name, contained in ${WD}
declare DATAEXPFLD=$(read_property "${CONFIGFILE}" "DATAEXPFLD") # The name of the experimental folder that we intend to analyze by running the script.
declare DATAEXPFILE=$(read_property "${CONFIGFILE}" "DATAEXPFILE") # The name of the target gene file that we try to analyze
declare DATAREFFILE=$(read_property "${CONFIGFILE}" "DATAREFFILE") # The name of the reference gene (housekeeping) file that we try to analyze
declare EXPLEVELS=""  # experimental categories
declare CTRLEVELS=""  # control categories
declare MONTECARLO_NPOINTS=$(read_property "${CONFIGFILE}" "MONTECARLO_NPOINTS") # number of points in the montecarlo simulation of the t-student variable difference.
declare QRANKS=$(read_property "${CONFIGFILE}" "QRANKS") # quantile ranks
declare WGRAP=$(read_property "${CONFIGFILE}" "WGRAP") # Width of the graph
declare WH=$(read_property "${CONFIGFILE}" "WH") # Width and height of the graph in inches separated by commas
declare LVE=""  # relative expression levels of the target gene (categories)
declare CVE=""  # relative expression levels of the reference gene (categories)
declare DIG=$(read_property "${CONFIGFILE}" "DIG") # quantile ranks



# PASSING ARGUMENTS --------------------
declare -ra ARGS=("$@")
for ((i=0; $i<${#ARGS[*]}; i++)); do
        if [ "${ARGS[${i}]}" == "--wd" ]; then
                WD="${ARGS[$((i+1))]}";
        elif [ "${ARGS[${i}]}" == "--analysis.folder" ]; then
                ANFLD="${ARGS[$((i+1))]}";
        elif [ "${ARGS[${i}]}" == "--data.folder" ]; then
                DATAFLD="${ARGS[$((i+1))]}";
        elif [ "${ARGS[${i}]}" == "--exp.folder" ]; then
		            DATAEXPFLD="${ARGS[$((i+1))]}";
	      elif [ "${ARGS[${i}]}" == "--exp.file" ]; then
		            DATAEXPFILE="${ARGS[$((i+1))]}";
	      elif [ "${ARGS[${i}]}" == "--ref.file" ]; then
		            DATAREFFILE="${ARGS[$((i+1))]}";
        elif [ "${ARGS[${i}]}" == "--exp.levels" ]; then
		            EXPLEVELS="${ARGS[$((i+1))]}";
        elif [ "${ARGS[${i}]}" == "--ctr.levels" ]; then
		            CTRLEVELS="${ARGS[$((i+1))]}";		
        elif [ "${ARGS[${i}]}" == "--montecarlo.n" ]; then
		            MONTECARLO_NPOINTS="${ARGS[$((i+1))]}";		
        elif [ "${ARGS[${i}]}" == "--qranks" ]; then
		            QRANKS="${ARGS[$((i+1))]}";		
		    elif [ "${ARGS[${i}]}" == "--wgraph" ]; then
		            WGRAP="${ARGS[$((i+1))]}";
	      elif [ "${ARGS[${i}]}" == "--wh" ]; then
		            WH="${ARGS[$((i+1))]}";
	      elif [ "${ARGS[${i}]}" == "--lve" ]; then
		            LVE="${ARGS[$((i+1))]}";		
	      elif [ "${ARGS[${i}]}" == "--cve" ]; then
		            CVE="${ARGS[$((i+1))]}";	
	      elif [ "${ARGS[${i}]}" == "--dig" ]; then
		            CVE="${ARGS[$((i+1))]}";				            
        fi                
done

DESTFLD="${WD}/${ANFLD}/${DATAEXPFLD}"
if [ ! -d "${DESTFLD}" ]; then
	echo "${BLUE_FORMAT}[INFO]${FORMAT_RESET} Creating analysis folder at ${DESTFLD} ..."
        mkdir -p "${DESTFLD}"
fi


R --vanilla --slave --file="${WD}/${ANFLD}/make_qpcr_analysis2.R" --args "${WD}" "${DATAFLD}" "${ANFLD}" "${DATAEXPFLD}"  "${DATAEXPFILE}" "${DATAREFFILE}" "${EXPLEVELS}" "${CTRLEVELS}" "${MONTECARLO_NPOINTS}" "${QRANKS}" "${WH}" "${WGRAP}"  "${LVE}" "${CVE}" "${DIG}"

exit 0
