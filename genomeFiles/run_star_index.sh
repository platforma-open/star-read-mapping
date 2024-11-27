#!/bin/bash

set -o nounset
set -o errexit

# Check for correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <Species> <Path_to_JSON_file_with_parameters> <Output_folder>"
    exit 1
fi

# Assign arguments to variables
: "${STAR_CMD:=STAR}"
SPECIES="${1}"
PARAMS_JSON_FILE="${2}"
OUTPUT_FOLDER="${3}"

echo "Running $SPECIES $PARAMS_JSON_FILE $OUTPUT_FOLDER "

# Automatically download STAR on Linux and Mac OS X
if ! command -v "${STAR_CMD}" >/dev/null; then
    _os="$(uname)"
    if [ "${_os}" != "Linux" ] && [ "${_os}" != "Darwin" ]; then
        echo "Executable '${STAR_CMD}' not found and this is not a Linux or Mac OS X OS."
        echo "Please, download and install STAR manually, add it to PATH or specify STAR_CMD env variable to point to STAR binary."
        exit 1
    fi

    _star_version="2.7.11b"

    if [ "${_os}" = "Darwin" ]; then
        STAR_CMD="downloads/STAR_${_star_version}/MacOSX_x86_64/STAR"
    else
        STAR_CMD="downloads/STAR_${_star_version}/Linux_x86_64_static/STAR"
    fi

    if ! [ -x "${STAR_CMD}" ]; then
        echo "No STAR executable found on host. Downloading '${_star_version}'..."
        (
            set -x
            mkdir -p "downloads"
            cd "downloads"
            wget --output-document "STAR.zip" \
                "https://github.com/alexdobin/STAR/releases/download/${_star_version}/STAR_${_star_version}.zip"
            unzip -o STAR.zip
        )
    fi

    STAR_CMD="$(realpath "${STAR_CMD}")"
fi

# Read species-specific parameters from JSON file
ASSEMBLY_VERSION=$(jq -r --arg species "$SPECIES" '.[$species].assembly_version' ${PARAMS_JSON_FILE})
READ_LENGTH=$(jq -r --arg species "$SPECIES" '.[$species].read_length' ${PARAMS_JSON_FILE})
GENOME_SA_INDEX_NBASES=$(jq -r --arg species "$SPECIES" '.[$species].genomeSAindexNbases' ${PARAMS_JSON_FILE})

# Calculate sjdbOverhang as read_length - 1
SJDB_OVERHANG=$((READ_LENGTH - 1))

# Get number of available CPU cores for threading
NUM_THREADS=$(nproc)

# Read genome and GTF URLs from the JSON file for the specified species
GENOME_URL=$(jq -r --arg species "$SPECIES" '.[$species].genome_url' ${PARAMS_JSON_FILE})
GTF_URL=$(jq -r --arg species "$SPECIES" '.[$species].gtf_url' ${PARAMS_JSON_FILE})

LOCAL_GENOME_FILE=$(jq -r --arg species "$SPECIES" '.[$species].local_genome_file' ${PARAMS_JSON_FILE})
LOCAL_GTF_FILE=$(jq -r --arg species "$SPECIES" '.[$species].local_gtf_file' ${PARAMS_JSON_FILE})

# We always should have either urls or local files
if [[ "$LOCAL_GENOME_FILE" != "null" && "$LOCAL_GTF_FILE" != "null" ]];then
  echo "Found local files: ${SPECIES}"
  LOCAL_FILES="true"
  REMOTE_FILES=""
elif [[ "$GENOME_URL" != "null" && "$GTF_URL" != "null" ]];then
  echo "Found url files: ${SPECIES}"
  REMOTE_FILES="true"
  LOCAL_FILES=""
else
 echo "Error: Neither GENOME/GTF Urls nor GENOME/GTF Local files are set."
 exit 1
fi

# Check if all required parameters were found
if [[ -z "$ASSEMBLY_VERSION" || -z "$READ_LENGTH" ]]; then
    echo "Error: Parameters for species '$SPECIES' not found in ${PARAMS_JSON_FILE}"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_FOLDER}
cd ${OUTPUT_FOLDER}


if [[ "$REMOTE_FILES" == "true" ]];then
  # Download and decompress genome DNA fasta file
  GENOME_FILENAME="genome.fa"
  echo "Downloading and decompressing genome for ${SPECIES} ${ASSEMBLY_VERSION}..."
  wget -O genome.fa.gz ${GENOME_URL}
  gunzip --force genome.fa.gz

  # Download genome annotation GTF file
  GTF_FILENAME="${SPECIES}_${ASSEMBLY_VERSION}_annotations.gtf"
  echo "Downloading and decompressing annotations for ${SPECIES} ${ASSEMBLY_VERSION}..."
  wget -O ${GTF_FILENAME}.gz ${GTF_URL}
  gunzip --force ${GTF_FILENAME}.gz
elif [[ "$LOCAL_FILES" == "true" ]];then
  CURRENT_DIR="${PARENT_SCRIPT_PATH}"
  GENOME_FILENAME="${CURRENT_DIR}/${LOCAL_GENOME_FILE}"
  GTF_FILENAME="${CURRENT_DIR}/${LOCAL_GTF_FILE}"
fi

# Generate Genome Index
echo "Generating genome index with sjdbOverhang=${SJDB_OVERHANG} and genomeSAindexNbases=${GENOME_SA_INDEX_NBASES}..."
"${STAR_CMD}" --runThreadN ${NUM_THREADS} --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ${GENOME_FILENAME} \
     --sjdbGTFfile ${GTF_FILENAME} --sjdbOverhang ${SJDB_OVERHANG} --genomeSAindexNbases ${GENOME_SA_INDEX_NBASES}

# Cleanup: Remove intermediate files
if [[ -f genome.fa ]];then
  rm genome.fa
fi
rm -rf _STARtmp/
# Rename GTF to standard name
mv "${GTF_FILENAME}" annotations.gtf

# echo "Output files: ${GTF_FILENAME} and ${TAR_FILENAME}"
echo "Genome index generation and packaging complete."
