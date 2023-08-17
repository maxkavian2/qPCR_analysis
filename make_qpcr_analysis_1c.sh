
#!/bin/bash

#
# Extracts and represent the data from QuantumStudio 3. 
# 
# This is the ordinary procedure. Per-sample technical replicas are reduced to their mean and
# the result is used for homogeneity tests (t-test and Kolmogorov-Smirnov test). 
# 
# 
# 
# USAGE: bash make_qpcr_analysis.sh --wd ${HOME} ...
# --analysis.folder qPCR_analysis
#	--data.folder qPCR_data
#	--exp.folder default
#	--exp.file default.csv
# --ref.file defaulf_ref.csv
#	--plate.ranks a1:b4,h8:k9,...
#	--categs standard,a1b4,h8k9,...
#	--categs.names "Standard,control, experimental ..."
#	--alpha 0.05
#	--wgraph 0.75
#	--wh 2,6
#	--delete.outliers TRUE
#	--delete.highsd FALSE
# 
# ARGUMENTS:
# --wd <string> : This is the top working directory (absolute path) were global folders for data and results can be found. 
# --analysis.folder <string> : name of the analysis folder where results will be stored. Avoid spaces for compatibility!
# --data.folder <string> : name of the data folder. Compatible files are in CSV format, with tab as separator and quoted text. Decimal character is a dot. Avoid spacies for compatibility!
# --exp.folder <string>: name of the experimental folder. Avoid spaces ...
# --exp.file <string>: name of the files from QS3. Avoid spaces ...
# --ref.file <string>: name of the files used as reference gene (e.g. housekeeping gene such as actin) for calculation purposes.
# --plate.ranks <string>: string containing the ranks of different plate categories separated by commas, in excel format.
# --categs <string>: categories associated to the plate ranks for internal data structure (defaults are give as excel format ranks without the colon character)
# --categs.names <string>: actual name of the categories, expressing conditions, genotypes, temperature, etc. meant to be displayed in the graphical output (chart)
# --alpha <float> : signification level for the CI.
# --wgraph <float> : width of each bar plot
# --wh <float,float> : width of the PDF file.
# --delete.outliers <boolean> : should outliers from QS3 files be removed?
# --delete.highsd <boolean> : should points with high SD be removed from QS3 files?
#
#
#  PDF representing the normalized values of the indicated categories as DeltaCT and Ratio, and several latex tables codes containing the results
#  of the t-tests and the KS-tests.
#
#  @AUTHOR: Maximo Sanchez-Aragon
#  @EMAIL: maxsanara@gmail.com
# 
# 



declare -r RED_FORMAT="$(printf '\033[31;1m')"
declare -r BLUE_FORMAT="$(printf '\033[34;1m')"
declare -r FORMAT_RESET="$(printf '\033[0m')"
declare -r CONFIGFILE_NAME="config_1c"
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
declare DATAEXPFILE=$(read_property "${CONFIGFILE}" "DATAEXPFILE") # The name of the experimental file that we try to analyze

# add line
declare DATAREFFILE=$(read_property "${CONFIGFILE}" "DATAREFFILE") # The name of the reference gene (housekeeping) file that we try to analyze

declare PLATERANKS="" # The plate ranks containing the different categories of data, in excel rank format, e.g. a1:b5, separated by commas.
declare CATEGS=$(echo "${PLATERANKS}" | tr -d ":" | tr -d " " | tr "[:lower:]" "[:upper:]") # Categories abbreviations used for plateranks. Default is a name that reminds the rank
declare CATEGS_NMS="${CATEGS}"  # Long names for the categories, by default they are the same as ${CATEGS}

declare STDCURVE_POS=$(read_property "${CONFIGFILE}" "STDCURVE_POS") # position of the standard curve
declare SKIPLINES=$(read_property "${CONFIGFILE}" "SKIPLINES") # number of lines that will be skipped while reading the quantum studio 3 file.
declare ALPHA=$(read_property "${CONFIGFILE}" "ALPHA") # Signification level for confidence intervals and statistical tests.
declare WGRAPH=$(read_property "${CONFIGFILE}" "WGRAPH") # Width of the graph
declare WH=$(read_property "${CONFIGFILE}" "WH") # Width and height of the graph in inches separated by commas
declare DELOUT=$(read_property "${CONFIGFILE}" "DELOUT") # Delete outliers before proceeding
declare DELHIGHSD=$(read_property "${CONFIGFILE}" "DELHIGHSD") # Delete high sd points before proceeding
declare EXCLUDING_FIELDS=$(read_property "${CONFIGFILE}" "EXCLUDING_FIELDS")  # Fields from the Quantum Studio 3 file excluded in the processing
declare FACTORING_FIELDS=$(read_property "${CONFIGFILE}" "FACTORING_FIELDS")  # Fields that will be factorized in R for compatibility




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
	elif [ "${ARGS[${i}]}" == "--plate.ranks" ]; then
		PLATERANKS="${ARGS[$((i+1))]}";
	elif [ "${ARGS[${i}]}" == "--categs" ]; then
		CATEGS="${ARGS[$((i+1))]}";
	elif [ "${ARGS[${i}]}" == "--categs.names" ]; then
		CATEGS_NMS="${ARGS[$((i+1))]}";	
	elif [ "${ARGS[${i}]}" == "--ref.file" ]; then
		            DATAREFFILE="${ARGS[$((i+1))]}";
	elif [ "${ARGS[${i}]}" == "--alpha" ]; then
		ALPHA="${ARGS[$((i+1))]}";
	elif [ "${ARGS[${i}]}" == "--wgraph" ]; then
		WGRAPH="${ARGS[$((i+1))]}";
	elif [ "${ARGS[${i}]}" == "--wh" ]; then
		ALPHA="${ARGS[$((i+1))]}";
	elif [ "${ARGS[${i}]}" == "--delete.outliers" ]; then
		DELOUT="${ARGS[$((i+1))]}";
	elif [ "${ARGS[${i}]}" == "--delete.highsd" ]; then
		DELHIGHSD="${ARGS[$((i+1))]}";
        fi                
done

DESTFLD="${WD}/${ANFLD}/${DATAEXPFLD}"
if [ ! -d "${DESTFLD}" ]; then
	echo "${BLUE_FORMAT}[INFO]${FORMAT_RESET} Creating analysis folder at ${DESTFLD} ..."
        mkdir -p "${DESTFLD}"
fi


R --vanilla --slave --file="${WD}/${ANFLD}/make_qpcr_analysis_1c.R" --args "${WD}" "${DATAFLD}" "${ANFLD}" "${DATAEXPFLD}" "${DATAEXPFILE}" "${DATAREFFILE}" "${PLATERANKS}" "${CATEGS}" "${CATEGS_NMS}" "${STDCURVE_POS}" "${SKIPLINES}" "${ALPHA}" "${WGRAPH}" "${WH}" "${DELOUT}" "${DELHIGHSD}" "${ALPHA}" "${EXCLUDING_FIELDS}" "${FACTORING_FIELDS}"



exit 0




