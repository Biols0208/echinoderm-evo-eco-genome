#!/usr/bin/env bash
set -e

# --- User and System Information (for logging) ---
# These are examples, the script will use the actual current user and time.
EXAMPLE_USER_LOGIN=$(whoami)
EXAMPLE_SCRIPT_EXECUTION_DATETIME=$(date +"%Y-%m-%d %H:%M:%S %Z")

CURRENT_USER_LOGIN_ACTUAL=$(whoami)
SCRIPT_EXECUTION_DATETIME_ACTUAL=$(date -u +"%Y-%m-%d %H:%M:%S UTC") # -u for UTC

echo "==================================================================================="
echo "SYNPHONI Multi-Node Pipeline Execution Log"
echo "Script Version: 2.5"
echo "User (Actual): $CURRENT_USER_LOGIN_ACTUAL"
echo "Execution Start Time (Actual): $SCRIPT_EXECUTION_DATETIME_ACTUAL"
echo "==================================================================================="
echo ""

# --- Python Executable and Script Directory Configuration ---
# or source activate conda env here
PYTHON_EXE="python"
SCRIPT_DIR="/soft/synphoni"
TOOLS_DIR="${SCRIPT_DIR}/tools"

echo "--- Configuration ---"
echo "Python Executable: $PYTHON_EXE"
echo "SYNPHONI Script Directory: $SCRIPT_DIR"
echo "SYNPHONI Tools Directory: $TOOLS_DIR"
echo ""

# --- Input File Configuration ---
CLUST_FILE="00_cluster/orthohmm_orthogroups.txt.2.tsv"
CHROM_FILES_PATTERN="00_cluster/chrom_files/*chrom"
SPECIES_TREE="00_cluster/species_tree.deal.nwk"

echo "Input Files:"
echo "  Orthology File (.clus): $CLUST_FILE"
echo "  Gene Location Files (.chrom pattern): $CHROM_FILES_PATTERN"
echo "  Species Tree File (.tre): $SPECIES_TREE"
echo ""

# --- Output Directory Configuration ---
CURRENT_RUN_BASE_OUTPUT_DIR="Synphoni_result_$(date +"%Y%m%d_%H%M%S")"
mkdir -p "$CURRENT_RUN_BASE_OUTPUT_DIR"
echo "Current Run Base Output Directory: $CURRENT_RUN_BASE_OUTPUT_DIR"
echo ""

# --- Node Selection for Analysis ---
#(((TRIAD:1,HOIHO:1)Placozoa:1,(ACRMI:1,(((LINAN:1,((ACHFU:1,POMCA:1)Apogastropoda:1,((PECMA:1,MIZYE:1)Pectinida:1,CRAGI:1)Pteriomorpha:1)Pleistomollusca:1)Lophotrochozoa:1,STRMA:1)Protostomia:1,(((ANNJA:1,((ACAPL:1,ASTRU:1)Echinozoa:1,STRPU:1):1)Echinodermata:1,(PTYFL:1,SACKO:1)Hemichordata:1)Ambulacraria:1,(BRAFL:1,(((HIPCO:1,MAYZE:1)Percomorpha:1,LEPOC:1)Neopterygii:1,((CANLUFA:1,HOMSA:1)Boreoeutheria:1,(PODMU:1,GALGA:1)Diapsida:1)Amniota:1)Osteichthyes:1)Chordata:1)Deuterostomia:1)Nephrozoa:1)Planulozoa:1)Parahoxozoa:1,EPHMU:1);
## example ALL_NODE_NAMES=("Nephrozoa" "Deuterostomia" "Ambulacraria" "Echinodermata" "Echinozoa")
ALL_NODE_NAMES=("budai" "Jipi" "haibaihe" "Asterozoa" "haishewei" "haixing" "Echinozoa" "Hemichordata" "haidan" "NODE18" "NODE19" "NODE20" "NODE14" "NODE16" "Haixing")
echo "Target Nodes for Analysis: ${ALL_NODE_NAMES[*]}"
if [ ${#ALL_NODE_NAMES[@]} -eq 0 ]; then
    echo "ERROR: ALL_NODE_NAMES array is empty. Please populate it with target node names."
    exit 1
fi
echo ""

# --- Resumption Configuration ---
# START_STEP: Define from which step to start/resume the analysis.
#             1: Start from the very beginning (Global Step 1).
#             2: Start from Node-specific Step 2.
#             2.5: Start from Node-specific Step 2.5 (nmax estimation).
#             3: Start from Node-specific Step 3.
#             4: Start from Node-specific Step 4.
#             5: Start from Node-specific Final Validation.
#             6: Start from Node-specific Optional Intervening Genes Analysis.
START_STEP="1" # <<< SET THIS TO THE DESIRED STARTING STEP

# PREVIOUS_BASE_OUTPUT_DIR: Path to the base output directory of a *previous* run.
PREVIOUS_BASE_OUTPUT_DIR="Synphoni_result"

echo "--- Resumption Configuration ---"
echo "Analysis will attempt to start from Step: $START_STEP"
if [ $(echo "$START_STEP > 1" | bc -l) -eq 1 ]; then
    if [ -z "$PREVIOUS_BASE_OUTPUT_DIR" ]; then
        echo "WARNING: START_STEP is > 1, but PREVIOUS_BASE_OUTPUT_DIR is not set."
        echo "         The script will expect Global Step 1 outputs (and any prior node-specific outputs if resuming > Step 2)"
        echo "         to be in the respective subdirectories of: '$CURRENT_RUN_BASE_OUTPUT_DIR'."
        EFFECTIVE_STEP1_OUTPUT_BASE_DIR="$CURRENT_RUN_BASE_OUTPUT_DIR"
        EFFECTIVE_NODE_PREREQ_BASE_DIR="$CURRENT_RUN_BASE_OUTPUT_DIR" # Node prereqs also from current run
    elif [ ! -d "$PREVIOUS_BASE_OUTPUT_DIR" ]; then
        echo "ERROR: PREVIOUS_BASE_OUTPUT_DIR '$PREVIOUS_BASE_OUTPUT_DIR' does not exist or is not a directory."
        echo "       Cannot resume from Step $START_STEP without valid prior Global Step 1 outputs."
        exit 1
    else
        echo "Resuming: Using Global Step 1 outputs (and potentially prior node-specific outputs) from PREVIOUS_BASE_OUTPUT_DIR: $PREVIOUS_BASE_OUTPUT_DIR"
        EFFECTIVE_STEP1_OUTPUT_BASE_DIR="$PREVIOUS_BASE_OUTPUT_DIR"
        EFFECTIVE_NODE_PREREQ_BASE_DIR="$PREVIOUS_BASE_OUTPUT_DIR" # Node prereqs also from previous run
    fi
else # START_STEP <= 1
    echo "Starting from Step 1 (or earlier). Global Step 1 outputs will be generated in: $CURRENT_RUN_BASE_OUTPUT_DIR/_GLOBAL_STEP1_OUTPUTS"
    EFFECTIVE_STEP1_OUTPUT_BASE_DIR="$CURRENT_RUN_BASE_OUTPUT_DIR"
    EFFECTIVE_NODE_PREREQ_BASE_DIR="$CURRENT_RUN_BASE_OUTPUT_DIR" # If starting fresh, node prereqs don't exist yet or will be made in current.
fi
echo ""

# --- Global Parameter Configuration for SYNPHONI Steps ---
# Step 1: step1_make_dists.py
#   -m, --maxpara: OGs with more paralogs from a single species than this will be excluded.
#                  Default: 100. Higher might include noisy, highly duplicated OGs.
MAX_PARA_STEP1=200

#   -s, --species_threshold: Min species an OG pair must be syntenic in to be retained.
#                            Default: 2. README suggests keeping this low for Step 1.
SPECIES_THRESHOLD_STEP1=2

# Step 2: step2_filter_pairs.py
#   -m, --species_threshold (referred to as 'm'): Min species in a phylogenetic clade
#                                                (ingroup, sister, outgroup sub-clades)
#                                                that must possess a syntenic OG pair for
#                                                that clade to be considered populated by it.
#                                                Default: 3. Higher 'm' is more stringent.
SPECIES_THRESHOLD_STEP2=3 # 'm' for step2

# Step 2.5: tools/step2.5_optimal_nmax.py (Optional nmax estimation)
# If true, attempts to estimate optimal nmax for each node.
# If false, or if estimation fails for a node, MANUAL_NMAX will be used.
RUN_STEP2_5=true

# MANUAL_NMAX: Distance threshold for Step 3 if nmax estimation is skipped or fails.
#              Default: 30 (from SYNPHONI README for step3). Example data uses 52.
#              This value is CRITICAL if not estimated.
MANUAL_NMAX=30 # Default if nmax estimation fails or is skipped in Step 2.5 .

# Step 4: step4_optimized.py
#   -l, --min_len: Min OG co-occurrences between two extant species blocks. Default: 3.
MIN_LEN_STEP4=3

#   -s, --min_shared: Min overlap coefficient between OG content of two extant blocks. Default: 0.5.
MIN_SHARED_STEP4="0.5"

#   -k, --clique_size: Min species in a multispecies block. Default: 3.
CLIQUE_SIZE_STEP4=3

#   -r, --min_community_coverage: Min percentage of ancestral OG community an extant block must possess. Default: 0.3.
MIN_COMMUNITY_COVERAGE_STEP4="0.3"

# chrom_clustering_method: k_clique; leiden
chrom_clustering_method="leiden"

cpu="10"
temp_path="$PWD/TMPdir"
[ -d $temp_path ] || mkdir $temp_path

BLOCKS_BY_NODE_SCRIPT="${TOOLS_DIR}/BlocksByNode.py" # <<< VERIFY THIS PATH

RUN_INTERVENING_GENES_ANALYSIS=true
echo "--- Global Parameters Set ---"
echo "MANUAL_NMAX (fallback for nmax estimation): $MANUAL_NMAX"
echo ""

# --- Define Paths for Global Step 1 Outputs using EFFECTIVE_STEP1_OUTPUT_BASE_DIR ---
STEP1_OUTPUT_DIR_LOCATION="${EFFECTIVE_STEP1_OUTPUT_BASE_DIR}/_GLOBAL_STEP1_OUTPUTS"
STEP1_PRIMARY_OUTPUT_PREFIX_FOR_SCRIPT="${STEP1_OUTPUT_DIR_LOCATION}/synphoni_run_global"
DIST_FILE="${STEP1_PRIMARY_OUTPUT_PREFIX_FOR_SCRIPT}/synphoni_run_global.dist"
CHROM_PICKLE="${STEP1_PRIMARY_OUTPUT_PREFIX_FOR_SCRIPT}/chromdata.pickle"

# --- Step 1: Build a global distance file (step1_make_dists.py) ---
if [ $(echo "$START_STEP <= 1" | bc -l) -eq 1 ]; then
    echo "==================================================================================="
    echo "Executing Step 1: Building Global Distance File (step1_make_dists.py)"
    echo "Outputs will be in: $STEP1_OUTPUT_DIR_LOCATION (within $CURRENT_RUN_BASE_OUTPUT_DIR as EFFECTIVE_STEP1_OUTPUT_BASE_DIR is $EFFECTIVE_STEP1_OUTPUT_BASE_DIR)"
    echo "==================================================================================="
    mkdir -p "$(dirname "$STEP1_PRIMARY_OUTPUT_PREFIX_FOR_SCRIPT")"

    CHROM_FILES_LIST=$(ls $CHROM_FILES_PATTERN 2>/dev/null)
    if [ -z "$CHROM_FILES_LIST" ]; then
        echo "ERROR: No .chrom files found matching pattern '$CHROM_FILES_PATTERN'. Cannot proceed with Step 1."
        exit 1
    fi

    echo "Running: $PYTHON_EXE ${SCRIPT_DIR}/step1_make_dists.py -c $CLUST_FILE $CHROM_FILES_LIST -o $STEP1_PRIMARY_OUTPUT_PREFIX_FOR_SCRIPT -m $MAX_PARA_STEP1 -s $SPECIES_THRESHOLD_STEP1"
    $PYTHON_EXE "${SCRIPT_DIR}/step1_make_dists.py" \
        -c "$CLUST_FILE" \
        $CHROM_FILES_LIST \
        -o "$STEP1_PRIMARY_OUTPUT_PREFIX_FOR_SCRIPT" \
        -m "$MAX_PARA_STEP1" \
        -s "$SPECIES_THRESHOLD_STEP1"


    if [ ! -f "$DIST_FILE" ]; then
        echo "ERROR: Step 1 FAILED. Global distance file '$DIST_FILE' was not created."
        exit 1
    fi
    if [ ! -f "$CHROM_PICKLE" ]; then
        echo "ERROR: Step 1 FAILED. Global pickle file '$CHROM_PICKLE' was not created."
        exit 1
    fi
    echo "Step 1 (Global) finished successfully."
    echo "  Global Distance File: $DIST_FILE"
    echo "  Global Chromosome Data Pickle: $CHROM_PICKLE"
    echo "==================================================================================="
    echo ""
else # START_STEP > 1
    echo "==================================================================================="
    echo "Skipping Step 1 (Global) execution as START_STEP = $START_STEP (> 1)."
    echo "Attempting to use existing Global Step 1 outputs from: $STEP1_OUTPUT_DIR_LOCATION"
    echo "  Expected Global Distance File: $DIST_FILE"
    echo "  Expected Global Chromosome Data Pickle: $CHROM_PICKLE"
    if [ ! -f "$DIST_FILE" ] || [ ! -f "$CHROM_PICKLE" ]; then
        echo "ERROR: Cannot skip Step 1. Required output file(s) from '$STEP1_OUTPUT_DIR_LOCATION' not found."
        exit 1
    fi
    echo "Required Global Step 1 files found. Proceeding to node-specific steps."
    echo "==================================================================================="
    echo ""
fi

# --- Function to attempt copying prerequisite files for a node ---
# Usage: copy_prerequisite_if_needed <source_file_path> <destination_file_path> <step_description> <node_name>
copy_prerequisite_if_needed() {
    local source_file="$1"
    local dest_file="$2"
    local step_description="$3"
    local node_name_for_log="$4"

    if [ ! -f "$dest_file" ]; then # If destination file doesn't exist in current run's node dir
        if [ -n "$PREVIOUS_BASE_OUTPUT_DIR" ] && [ -d "$PREVIOUS_BASE_OUTPUT_DIR" ] && [ -f "$source_file" ]; then
            echo "INFO for Node '$node_name_for_log': Prerequisite for $step_description ('$(basename "$dest_file")') not found in current run."
            echo "      Attempting to copy from PREVIOUS_BASE_OUTPUT_DIR: $source_file"
            mkdir -p "$(dirname "$dest_file")" # Ensure destination directory exists in current run
            cp "$source_file" "$dest_file"
            if [ $? -eq 0 ]; then
                echo "      Successfully copied to: $dest_file"
            else
                echo "      ERROR: Failed to copy prerequisite file '$source_file' for Node '$node_name_for_log'."
                return 1 # Indicate failure
            fi
        else
            echo "ERROR for Node '$node_name_for_log': Prerequisite for $step_description ('$(basename "$dest_file")') not found in current run."
            if [ -n "$PREVIOUS_BASE_OUTPUT_DIR" ] && [ -d "$PREVIOUS_BASE_OUTPUT_DIR" ]; then
                 echo "       And could not be found in PREVIOUS_BASE_OUTPUT_DIR (Source looked for: '$source_file' - Exists? $(test -f "$source_file" && echo Yes || echo No))."
            else
                 echo "       PREVIOUS_BASE_OUTPUT_DIR was not set or invalid, so no copy attempted from there."
            fi
            echo "       Please ensure this file exists in the expected location or adjust START_STEP to generate it."
            return 1 # Indicate failure
        fi
    else
        echo "INFO for Node '$node_name_for_log': Prerequisite for $step_description ('$(basename "$dest_file")') already exists in current run: $dest_file"
    fi
    return 0 # Indicate success or file already present
}


# --- Loop through each specified node for Steps 2 onwards ---
echo "Starting analysis loop for each specified node (from step $START_STEP)..."
for CURRENT_NODE_NAME in "${ALL_NODE_NAMES[@]}"; do
    echo ""
    echo "###################################################################################"
    echo "Processing Node: $CURRENT_NODE_NAME (Attempting from Step $START_STEP)"
    echo "###################################################################################"

    NODE_SPECIFIC_OUTPUT_DIR_CURRENT_RUN="${CURRENT_RUN_BASE_OUTPUT_DIR}/${CURRENT_NODE_NAME}"
    mkdir -p "$NODE_SPECIFIC_OUTPUT_DIR_CURRENT_RUN"
    echo "Node-specific output directory for this run: $NODE_SPECIFIC_OUTPUT_DIR_CURRENT_RUN"

    # Define paths for expected outputs IN THE CURRENT RUN'S DIRECTORY
    BASENAME_DIST_FILE_GLOBAL=$(basename "$DIST_FILE" .dist)
    NODE_DIST_FILE_FOR_STEP3_CURRENT_RUN="${NODE_SPECIFIC_OUTPUT_DIR_CURRENT_RUN}/${BASENAME_DIST_FILE_GLOBAL}.m_${SPECIES_THRESHOLD_STEP2}.${CURRENT_NODE_NAME}.dist"
    OG_COMMUS_OUTPUT_PREFIX_NODE_CURRENT_RUN="${NODE_SPECIFIC_OUTPUT_DIR_CURRENT_RUN}/${CURRENT_NODE_NAME}.og_commus"
    OG_COMMUS_CSV_NODE_CURRENT_RUN="${OG_COMMUS_OUTPUT_PREFIX_NODE_CURRENT_RUN}.csv"
    OG_COMMUS_GPICKLE_NODE_CURRENT_RUN="${OG_COMMUS_OUTPUT_PREFIX_NODE_CURRENT_RUN}.gpickle"
    BLOCKS_OUTPUT_PREFIX_NODE_CURRENT_RUN="${NODE_SPECIFIC_OUTPUT_DIR_CURRENT_RUN}/${CURRENT_NODE_NAME}_blocks"
    STEP4_CLUSTERS_FILE_NODE_CURRENT_RUN="${BLOCKS_OUTPUT_PREFIX_NODE_CURRENT_RUN}.len${MIN_LEN_STEP4}.ol${MIN_SHARED_STEP4}.clusters"
    STEP4_SYNT_FILE_NODE_CURRENT_RUN="${BLOCKS_OUTPUT_PREFIX_NODE_CURRENT_RUN}.len${MIN_LEN_STEP4}.ol${MIN_SHARED_STEP4}.synt"
    VALIDATED_SYNT_FILE_NODE_CURRENT_RUN="${NODE_SPECIFIC_OUTPUT_DIR_CURRENT_RUN}/${CURRENT_NODE_NAME}_blocks.len${MIN_LEN_STEP4}.ol${MIN_SHARED_STEP4}.taxonomy_filtered.synt"
    VALIDATED_CLUSTERS_FILE_NODE_CURRENT_RUN="${NODE_SPECIFIC_OUTPUT_DIR_CURRENT_RUN}/${CURRENT_NODE_NAME}_blocks.len${MIN_LEN_STEP4}.ol${MIN_SHARED_STEP4}.taxonomy_filtered.clusters"
    NMAX_INFLECTION_FILE_NODE_CURRENT_RUN="${NODE_SPECIFIC_OUTPUT_DIR_CURRENT_RUN}/nmax_${CURRENT_NODE_NAME}.inflection.csv"
    INTERVENING_GENES_OUTPUT_PREFIX_NODE_CURRENT_RUN="${NODE_SPECIFIC_OUTPUT_DIR_CURRENT_RUN}/${CURRENT_NODE_NAME}_intervening_genes"
    INTERVENING_GENES_TSV_NODE_CURRENT_RUN="${INTERVENING_GENES_OUTPUT_PREFIX_NODE_CURRENT_RUN}.tsv"


    # Define paths for prerequisite files FROM THE PREVIOUS RUN'S DIRECTORY (using EFFECTIVE_NODE_PREREQ_BASE_DIR)
    NODE_SPECIFIC_PREREQ_BASE_DIR="${EFFECTIVE_NODE_PREREQ_BASE_DIR}/${CURRENT_NODE_NAME}" # Base for this node's prereqs

    NODE_DIST_FILE_FOR_STEP3_PREREQ="${NODE_SPECIFIC_PREREQ_BASE_DIR}/${BASENAME_DIST_FILE_GLOBAL}.m_${SPECIES_THRESHOLD_STEP2}.${CURRENT_NODE_NAME}.dist"
    OG_COMMUS_CSV_NODE_PREREQ="${NODE_SPECIFIC_PREREQ_BASE_DIR}/${CURRENT_NODE_NAME}.og_commus.csv"
    OG_COMMUS_GPICKLE_NODE_PREREQ="${NODE_SPECIFIC_PREREQ_BASE_DIR}/${CURRENT_NODE_NAME}.og_commus.gpickle"
    STEP4_CLUSTERS_FILE_NODE_PREREQ="${NODE_SPECIFIC_PREREQ_BASE_DIR}/${CURRENT_NODE_NAME}_blocks.len${MIN_LEN_STEP4}.ol${MIN_SHARED_STEP4}.clusters"
    STEP4_SYNT_FILE_NODE_PREREQ="${NODE_SPECIFIC_PREREQ_BASE_DIR}/${CURRENT_NODE_NAME}_blocks.len${MIN_LEN_STEP4}.ol${MIN_SHARED_STEP4}.synt"
    VALIDATED_SYNT_FILE_NODE_PREREQ="${NODE_SPECIFIC_PREREQ_BASE_DIR}/${CURRENT_NODE_NAME}_blocks.len${MIN_LEN_STEP4}.ol${MIN_SHARED_STEP4}.taxonomy_filtered.synt"
    VALIDATED_CLUSTERS_FILE_NODE_PREREQ="${NODE_SPECIFIC_PREREQ_BASE_DIR}/${CURRENT_NODE_NAME}_blocks.len${MIN_LEN_STEP4}.ol${MIN_SHARED_STEP4}.taxonomy_filtered.clusters"
    NMAX_INFLECTION_FILE_NODE_PREREQ="${NODE_SPECIFIC_PREREQ_BASE_DIR}/nmax_${CURRENT_NODE_NAME}.inflection.csv"


    # --- Step 2: Infer ancestral OG networks ---
    if [ $(echo "$START_STEP <= 2" | bc -l) -eq 1 ]; then
        echo ""
        echo "-----------------------------------------------------------------------------------"
        echo "Executing Step 2 for Node '$CURRENT_NODE_NAME': Inferring Ancestral OG Network"
        echo "Output will be in: $NODE_SPECIFIC_OUTPUT_DIR_CURRENT_RUN"
        echo "-----------------------------------------------------------------------------------"
        TEMP_NODE_DIST_FILE_IN_CWD="${BASENAME_DIST_FILE_GLOBAL}.m_${SPECIES_THRESHOLD_STEP2}.${CURRENT_NODE_NAME}.dist"
        echo "Running: $PYTHON_EXE ${SCRIPT_DIR}/step2_filter_pairs.py -d $DIST_FILE -s $SPECIES_TREE -n $CURRENT_NODE_NAME -m $SPECIES_THRESHOLD_STEP2"
        $PYTHON_EXE "${SCRIPT_DIR}/step2_filter_pairs.py" -d "$DIST_FILE" -s "$SPECIES_TREE" -n "$CURRENT_NODE_NAME" -m "$SPECIES_THRESHOLD_STEP2"
        if [ -f "$TEMP_NODE_DIST_FILE_IN_CWD" ]; then
            mv "$TEMP_NODE_DIST_FILE_IN_CWD" "$NODE_DIST_FILE_FOR_STEP3_CURRENT_RUN"
        elif [ ! -f "$NODE_DIST_FILE_FOR_STEP3_CURRENT_RUN" ]; then # Check if already moved or created there
            echo "ERROR: Step 2 FAILED for Node '$CURRENT_NODE_NAME'. Output file not found." ; continue ; fi
        echo "Step 2 for Node '$CURRENT_NODE_NAME' finished."
    else # START_STEP > 2
        echo "Skipping Step 2 execution for Node '$CURRENT_NODE_NAME' (START_STEP=$START_STEP)."
        copy_prerequisite_if_needed "$NODE_DIST_FILE_FOR_STEP3_PREREQ" "$NODE_DIST_FILE_FOR_STEP3_CURRENT_RUN" "Step 2 output (for Step 2.5/3 input)" "$CURRENT_NODE_NAME"
        if [ $? -ne 0 ] || [ ! -f "$NODE_DIST_FILE_FOR_STEP3_CURRENT_RUN" ]; then
             echo "ERROR: Prerequisite for Step 2.5/3 (from Step 2) missing for '$CURRENT_NODE_NAME'. Skipping node."; continue;
        fi
    fi
    echo "  Using Step 2 output: $NODE_DIST_FILE_FOR_STEP3_CURRENT_RUN"

    # --- Step 2.5: Estimate optimal nmax ---
    NMAX_VALUE=$MANUAL_NMAX # Reset for each node before potential estimation
    if [ $(echo "$START_STEP <= 2.5" | bc -l) -eq 1 ]; then
        if [ "$RUN_STEP2_5" = true ]; then
            echo ""
            echo "-----------------------------------------------------------------------------------"
            echo "Executing Step 2.5 for Node '$CURRENT_NODE_NAME': Estimating Optimal nmax"
            echo "Using input: $NODE_DIST_FILE_FOR_STEP3_CURRENT_RUN"
            echo "Output will be in: $NODE_SPECIFIC_OUTPUT_DIR_CURRENT_RUN (nmax_...)"
            echo "-----------------------------------------------------------------------------------"
            NMAX_OPTIMAL_OUTPUT_PREFIX_NODE_CURRENT_RUN="${NODE_SPECIFIC_OUTPUT_DIR_CURRENT_RUN}/nmax_${CURRENT_NODE_NAME}"
            STEP2_5_PARENT_DIR_FOR_NODE="${NODE_SPECIFIC_OUTPUT_DIR_CURRENT_RUN}/_step2.5_temp_input" # Temp dir for script structure
            STEP2_5_ACTUAL_INPUT_DIR_FOR_SCRIPT="${STEP2_5_PARENT_DIR_FOR_NODE}/${CURRENT_NODE_NAME}/m_${SPECIES_THRESHOLD_STEP2}"
            mkdir -p "$STEP2_5_ACTUAL_INPUT_DIR_FOR_SCRIPT"
            cp "$NODE_DIST_FILE_FOR_STEP3_CURRENT_RUN" "$STEP2_5_ACTUAL_INPUT_DIR_FOR_SCRIPT/"

            echo "Running: $PYTHON_EXE ${TOOLS_DIR}/step2.5_optimal_nmax.py $STEP2_5_PARENT_DIR_FOR_NODE/ -o $NMAX_OPTIMAL_OUTPUT_PREFIX_NODE_CURRENT_RUN"
            $PYTHON_EXE "${TOOLS_DIR}/step2.5_optimal_nmax.py" "$STEP2_5_PARENT_DIR_FOR_NODE/" -o "$NMAX_OPTIMAL_OUTPUT_PREFIX_NODE_CURRENT_RUN"
            # rm -rf "$STEP2_5_PARENT_DIR_FOR_NODE" # Optional cleanup of temp input structure

            if [ -f "$NMAX_INFLECTION_FILE_NODE_CURRENT_RUN" ]; then
                RAW_PARSED_NMAX_VALUE=$(awk -F',' -v node="$CURRENT_NODE_NAME" -v m_val="$SPECIES_THRESHOLD_STEP2" '($1==node && $2==m_val){print $3;exit}' "$NMAX_INFLECTION_FILE_NODE_CURRENT_RUN" | sed 's/\.$//')
                CLEANED_PARSED_NMAX_VALUE=$(echo "$RAW_PARSED_NMAX_VALUE" | tr -d '\r' | sed 's/[^0-9.]*$//' | sed 's/\.$//')
                if [ -n "$CLEANED_PARSED_NMAX_VALUE" ] && [[ "$CLEANED_PARSED_NMAX_VALUE" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
                    NMAX_VALUE=$(printf "%.0f" "$CLEANED_PARSED_NMAX_VALUE")
                    echo "Step 2.5 estimated nmax for '$CURRENT_NODE_NAME': $NMAX_VALUE"
                else
                    echo "WARNING: Could not parse nmax from Step 2.5 output for '$CURRENT_NODE_NAME' (Raw: '$RAW_PARSED_NMAX_VALUE', Cleaned: '$CLEANED_PARSED_NMAX_VALUE'). Using manual: $MANUAL_NMAX"
                    NMAX_VALUE=$MANUAL_NMAX
                fi
            else
                echo "WARNING: Step 2.5 output file '$NMAX_INFLECTION_FILE_NODE_CURRENT_RUN' not found for '$CURRENT_NODE_NAME'. Using manual nmax: $MANUAL_NMAX"
                NMAX_VALUE=$MANUAL_NMAX
            fi
        else # RUN_STEP2_5 is false
            echo "Skipping Step 2.5 execution (RUN_STEP2_5=false) for '$CURRENT_NODE_NAME'. Using manual nmax: $MANUAL_NMAX"
        fi
    else # START_STEP > 2.5, try to use existing nmax file if RUN_STEP2_5 was true
        echo "Skipping Step 2.5 execution for Node '$CURRENT_NODE_NAME' (START_STEP=$START_STEP)."
        if [ "$RUN_STEP2_5" = true ]; then # Only try to read/copy if it would have been run
            copy_prerequisite_if_needed "$NMAX_INFLECTION_FILE_NODE_PREREQ" "$NMAX_INFLECTION_FILE_NODE_CURRENT_RUN" "Step 2.5 output (nmax inflection file)" "$CURRENT_NODE_NAME"
            # If copy_prerequisite_if_needed failed, $? will be non-zero, but we might still proceed with MANUAL_NMAX if file is not there after attempt.
            if [ -f "$NMAX_INFLECTION_FILE_NODE_CURRENT_RUN" ]; then
                RAW_PARSED_NMAX_VALUE=$(awk -F',' -v node="$CURRENT_NODE_NAME" -v m_val="$SPECIES_THRESHOLD_STEP2" '($1==node && $2==m_val){print $3;exit}' "$NMAX_INFLECTION_FILE_NODE_CURRENT_RUN" | sed 's/\.$//')
                CLEANED_PARSED_NMAX_VALUE=$(echo "$RAW_PARSED_NMAX_VALUE" | tr -d '\r' | sed 's/[^0-9.]*$//' | sed 's/\.$//')
                if [ -n "$CLEANED_PARSED_NMAX_VALUE" ] && [[ "$CLEANED_PARSED_NMAX_VALUE" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
                    NMAX_VALUE=$(printf "%.0f" "$CLEANED_PARSED_NMAX_VALUE")
                    echo "Using previously available/copied nmax for '$CURRENT_NODE_NAME': $NMAX_VALUE"
                else
                    echo "WARNING: Could not parse previously available/copied nmax for '$CURRENT_NODE_NAME' (Raw: '$RAW_PARSED_NMAX_VALUE', Cleaned: '$CLEANED_PARSED_NMAX_VALUE'). Using manual: $MANUAL_NMAX"
                    # NMAX_VALUE remains $MANUAL_NMAX
                fi
            else
                echo "WARNING: Previous Step 2.5 output file ('$NMAX_INFLECTION_FILE_NODE_CURRENT_RUN') not found or copied for '$CURRENT_NODE_NAME'. Using manual nmax: $MANUAL_NMAX"
                # NMAX_VALUE remains $MANUAL_NMAX
            fi
        else # RUN_STEP2_5 is false
             echo "Using manual nmax: $MANUAL_NMAX (RUN_STEP2_5 is false)."
        fi
    fi
    echo "  nmax value for Step 3 (Node '$CURRENT_NODE_NAME'): $NMAX_VALUE"

    # --- Step 3: Isolate OG sets ---
    if [ $(echo "$START_STEP <= 3" | bc -l) -eq 1 ]; then
        echo ""
        echo "-----------------------------------------------------------------------------------"
        echo "Executing Step 3 for Node '$CURRENT_NODE_NAME' (nmax=$NMAX_VALUE)..."
        echo "Input: $NODE_DIST_FILE_FOR_STEP3_CURRENT_RUN"
        echo "Output prefix: $OG_COMMUS_OUTPUT_PREFIX_NODE_CURRENT_RUN"
        echo "-----------------------------------------------------------------------------------"
        $PYTHON_EXE "${SCRIPT_DIR}/step3_find_og_commus.py" "$NODE_DIST_FILE_FOR_STEP3_CURRENT_RUN" -n "$NMAX_VALUE" -o "$OG_COMMUS_OUTPUT_PREFIX_NODE_CURRENT_RUN"
        if [ ! -f "$OG_COMMUS_CSV_NODE_CURRENT_RUN" ] || [ ! -f "$OG_COMMUS_GPICKLE_NODE_CURRENT_RUN" ]; then
            echo "ERROR: Step 3 FAILED for Node '$CURRENT_NODE_NAME'."; continue; fi
        echo "Step 3 for Node '$CURRENT_NODE_NAME' finished."
    else # START_STEP > 3
        echo "Skipping Step 3 execution for Node '$CURRENT_NODE_NAME' (START_STEP=$START_STEP)."
        copy_prerequisite_if_needed "$OG_COMMUS_CSV_NODE_PREREQ" "$OG_COMMUS_CSV_NODE_CURRENT_RUN" "Step 3 output (CSV)" "$CURRENT_NODE_NAME"
        local csv_copied_ok=$?
        copy_prerequisite_if_needed "$OG_COMMUS_GPICKLE_NODE_PREREQ" "$OG_COMMUS_GPICKLE_NODE_CURRENT_RUN" "Step 3 output (GPICKLE)" "$CURRENT_NODE_NAME"
        local gpickle_copied_ok=$?
        if [ $csv_copied_ok -ne 0 ] || [ $gpickle_copied_ok -ne 0 ] || [ ! -f "$OG_COMMUS_CSV_NODE_CURRENT_RUN" ] || [ ! -f "$OG_COMMUS_GPICKLE_NODE_CURRENT_RUN" ]; then
            echo "ERROR: Prerequisites for Step 4 (from Step 3) missing for '$CURRENT_NODE_NAME'. Skipping node."; continue; fi
    fi
    echo "  Using Step 3 outputs: $OG_COMMUS_CSV_NODE_CURRENT_RUN, $OG_COMMUS_GPICKLE_NODE_CURRENT_RUN"

    # --- Step 4: Recover ancestral microsyntenic blocks ---
    if [ $(echo "$START_STEP <= 4" | bc -l) -eq 1 ]; then
        echo ""
        echo "-----------------------------------------------------------------------------------"
        echo "Executing Step 4 for Node '$CURRENT_NODE_NAME'..."
        echo "Inputs: $OG_COMMUS_CSV_NODE_CURRENT_RUN, $OG_COMMUS_GPICKLE_NODE_CURRENT_RUN, $CHROM_PICKLE"
        echo "Output prefix: $BLOCKS_OUTPUT_PREFIX_NODE_CURRENT_RUN"
        echo "-----------------------------------------------------------------------------------"
        path_new=$PWD
        cd "${SCRIPT_DIR}"
        $PYTHON_EXE -m "memory_profiler" "step4_optimized.py" -g "$path_new/$OG_COMMUS_GPICKLE_NODE_CURRENT_RUN" -c "$path_new/$CHROM_PICKLE" -o "$path_new/$BLOCKS_OUTPUT_PREFIX_NODE_CURRENT_RUN" -l $MIN_LEN_STEP4 -s $MIN_SHARED_STEP4 -k $CLIQUE_SIZE_STEP4 -r $MIN_COMMUNITY_COVERAGE_STEP4 --num-workers $cpu --temp-path "$temp_path" -m "$chrom_clustering_method" "$path_new/$OG_COMMUS_CSV_NODE_CURRENT_RUN"
        cd "${path_new}"
        if [ ! -f "$STEP4_CLUSTERS_FILE_NODE_CURRENT_RUN" ] || [ ! -f "$STEP4_SYNT_FILE_NODE_CURRENT_RUN" ]; then
            echo "ERROR: Step 4 FAILED for Node '$CURRENT_NODE_NAME'."; continue; fi
        echo "Step 4 for Node '$CURRENT_NODE_NAME' finished."
    else # START_STEP > 4
        echo "Skipping Step 4 execution for Node '$CURRENT_NODE_NAME' (START_STEP=$START_STEP)."
        copy_prerequisite_if_needed "$STEP4_CLUSTERS_FILE_NODE_PREREQ" "$STEP4_CLUSTERS_FILE_NODE_CURRENT_RUN" "Step 4 output (Clusters)" "$CURRENT_NODE_NAME"
        local clusters_copied_ok=$?
        copy_prerequisite_if_needed "$STEP4_SYNT_FILE_NODE_PREREQ" "$STEP4_SYNT_FILE_NODE_CURRENT_RUN" "Step 4 output (Synt)" "$CURRENT_NODE_NAME"
        local synt_copied_ok=$?
        if [ $clusters_copied_ok -ne 0 ] || [ $synt_copied_ok -ne 0 ] || [ ! -f "$STEP4_CLUSTERS_FILE_NODE_CURRENT_RUN" ] || [ ! -f "$STEP4_SYNT_FILE_NODE_CURRENT_RUN" ]; then
            echo "ERROR: Prerequisites for Step 5 (from Step 4) missing for '$CURRENT_NODE_NAME'. Skipping node."; continue; fi
    fi
    echo "  Using Step 4 outputs: $STEP4_CLUSTERS_FILE_NODE_CURRENT_RUN, $STEP4_SYNT_FILE_NODE_CURRENT_RUN"

    # --- Step 5: Final Validation ---
    if [ $(echo "$START_STEP <= 5" | bc -l) -eq 1 ]; then
        echo ""
        echo "-----------------------------------------------------------------------------------"
        echo "Executing Step 5 (Final Validation) for Node '$CURRENT_NODE_NAME' (BlocksByNode.py)"
        echo "Inputs: $STEP4_CLUSTERS_FILE_NODE_CURRENT_RUN, $STEP4_SYNT_FILE_NODE_CURRENT_RUN"
        echo "Outputs: $VALIDATED_SYNT_FILE_NODE_CURRENT_RUN, $VALIDATED_CLUSTERS_FILE_NODE_CURRENT_RUN"
        echo "-----------------------------------------------------------------------------------"
        if [ ! -f "$BLOCKS_BY_NODE_SCRIPT" ]; then echo "ERROR: BlocksByNode.py not found at '$BLOCKS_BY_NODE_SCRIPT'. Skipping validation.";
        else
            echo "Running validation for .synt file..."
            dos2unix "$STEP4_CLUSTERS_FILE_NODE_CURRENT_RUN"
            dos2unix "$STEP4_SYNT_FILE_NODE_CURRENT_RUN"
            cp "$STEP4_CLUSTERS_FILE_NODE_CURRENT_RUN" "${STEP4_CLUSTERS_FILE_NODE_CURRENT_RUN}.bak"
            cut -f 1,3- "${STEP4_CLUSTERS_FILE_NODE_CURRENT_RUN}.bak" > "$STEP4_CLUSTERS_FILE_NODE_CURRENT_RUN"

##  -r {short,clusters_list,blocks_list,tree_ASCII,tree_NH}, --report {short,clusters_list,blocks_list,tree_ASCII,tree_NH}
##                      Type of report printed to STDOUT [Default: short] : - "short": number of blocks per node. - "clusters_list": filters one multi-species block per
##                      line, cluster IDs (field 1), ancestral/novel nodes of the block (field 2), species list (field 3), block_ids (field 4+). For getting a *.clusters
##                      file with only a subset, pipe SyntByNode output to `cut -f1,4-`. - "blocks_list": blocks within the filtered multi-species blocks - "tree_ASCII":
##                      block count per specified node on an ASCII tree - "tree_NH": block count per specified node on a newick tree

            $PYTHON_EXE "$BLOCKS_BY_NODE_SCRIPT" -c "$STEP4_CLUSTERS_FILE_NODE_CURRENT_RUN" -b "$STEP4_SYNT_FILE_NODE_CURRENT_RUN" -s "$SPECIES_TREE" -n "$CURRENT_NODE_NAME" -m "$SPECIES_THRESHOLD_STEP2" -r short -t total | cut -f2- > "${VALIDATED_SYNT_FILE_NODE_CURRENT_RUN}.stat"
            $PYTHON_EXE "$BLOCKS_BY_NODE_SCRIPT" -c "$STEP4_CLUSTERS_FILE_NODE_CURRENT_RUN" -b "$STEP4_SYNT_FILE_NODE_CURRENT_RUN" -s "$SPECIES_TREE" -n "$CURRENT_NODE_NAME" -m "$SPECIES_THRESHOLD_STEP2" -r blocks_list -t total | cut -f2- > "$VALIDATED_SYNT_FILE_NODE_CURRENT_RUN"
            echo "Running validation for .clusters file..."
            $PYTHON_EXE "$BLOCKS_BY_NODE_SCRIPT" -c "$STEP4_CLUSTERS_FILE_NODE_CURRENT_RUN" -b "$STEP4_SYNT_FILE_NODE_CURRENT_RUN" -s "$SPECIES_TREE" -n "$CURRENT_NODE_NAME" -m "$SPECIES_THRESHOLD_STEP2" -r clusters_list -t total | cut -f1,3- > "$VALIDATED_CLUSTERS_FILE_NODE_CURRENT_RUN"

            # Check if validated files are non-empty (empty might be a valid result)
            if [ -f "$VALIDATED_SYNT_FILE_NODE_CURRENT_RUN" ] && [ ! -s "$VALIDATED_SYNT_FILE_NODE_CURRENT_RUN" ]; then echo "WARNING: Validated .synt file for '$CURRENT_NODE_NAME' is empty."; fi
            if [ -f "$VALIDATED_CLUSTERS_FILE_NODE_CURRENT_RUN" ] && [ ! -s "$VALIDATED_CLUSTERS_FILE_NODE_CURRENT_RUN" ]; then echo "WARNING: Validated .clusters file for '$CURRENT_NODE_NAME' is empty."; fi
            echo "Step 5 (Final Validation) for Node '$CURRENT_NODE_NAME' finished."
        fi
    else # START_STEP > 5
        echo "Skipping Step 5 execution for Node '$CURRENT_NODE_NAME' (START_STEP=$START_STEP)."
        copy_prerequisite_if_needed "$VALIDATED_SYNT_FILE_NODE_PREREQ" "$VALIDATED_SYNT_FILE_NODE_CURRENT_RUN" "Step 5 output (Validated Synt)" "$CURRENT_NODE_NAME"
        copy_prerequisite_if_needed "$VALIDATED_CLUSTERS_FILE_NODE_PREREQ" "$VALIDATED_CLUSTERS_FILE_NODE_CURRENT_RUN" "Step 5 output (Validated Clusters)" "$CURRENT_NODE_NAME"
        # No critical error if these are missing for step 6, as step 6 has fallbacks and is optional.
        if [ "$RUN_INTERVENING_GENES_ANALYSIS" = true ]; then
             if [ ! -f "$VALIDATED_SYNT_FILE_NODE_CURRENT_RUN" ] || [ ! -f "$VALIDATED_CLUSTERS_FILE_NODE_CURRENT_RUN" ]; then
                 echo "Note for Step 6: Validated files from Step 5 (or previous run) are expected for intervening gene analysis."
             fi
        fi
    fi
    echo "  Using Step 5 outputs (Validated): $VALIDATED_SYNT_FILE_NODE_CURRENT_RUN, $VALIDATED_CLUSTERS_FILE_NODE_CURRENT_RUN"

    # --- Step 6: Optional Intervening Genes Analysis ---
    if [ $(echo "$START_STEP <= 6" | bc -l) -eq 1 ]; then
        if [ "$RUN_INTERVENING_GENES_ANALYSIS" = true ]; then
            echo ""
            echo "-----------------------------------------------------------------------------------"
            echo "Executing Step 6 (Optional Intervening Genes Analysis) for Node '$CURRENT_NODE_NAME'"
            echo "Output prefix: $INTERVENING_GENES_OUTPUT_PREFIX_NODE_CURRENT_RUN"
            echo "-----------------------------------------------------------------------------------"
            SYNT_FOR_ANALYSIS=$VALIDATED_SYNT_FILE_NODE_CURRENT_RUN
            CLUSTERS_FOR_ANALYSIS=$VALIDATED_CLUSTERS_FILE_NODE_CURRENT_RUN

            # Fallback logic if validated files are missing or empty (use Step 4 outputs instead)
            if [ ! -f "$SYNT_FOR_ANALYSIS" ] || [ ! -s "$SYNT_FOR_ANALYSIS" ]; then
                echo "INFO for Step 6: Validated .synt file ('$(basename "$SYNT_FOR_ANALYSIS")') is empty or missing. Using Step 4 .synt file ('$(basename "$STEP4_SYNT_FILE_NODE_CURRENT_RUN")') instead."
                SYNT_FOR_ANALYSIS=$STEP4_SYNT_FILE_NODE_CURRENT_RUN
            fi
            if [ ! -f "$CLUSTERS_FOR_ANALYSIS" ] || [ ! -s "$CLUSTERS_FOR_ANALYSIS" ]; then
                echo "INFO for Step 6: Validated .clusters file ('$(basename "$CLUSTERS_FOR_ANALYSIS")') is empty or missing. Using Step 4 .clusters file ('$(basename "$STEP4_CLUSTERS_FILE_NODE_CURRENT_RUN")') instead."
                CLUSTERS_FOR_ANALYSIS=$STEP4_CLUSTERS_FILE_NODE_CURRENT_RUN
            fi

            if [ -f "$SYNT_FOR_ANALYSIS" ] && [ -f "$CLUSTERS_FOR_ANALYSIS" ]; then
                echo "Running: $PYTHON_EXE ${SCRIPT_DIR}/analysis_intervening_genes.py -c $CHROM_PICKLE -sy $SYNT_FOR_ANALYSIS -ms $CLUSTERS_FOR_ANALYSIS -o $INTERVENING_GENES_OUTPUT_PREFIX_NODE_CURRENT_RUN"
                $PYTHON_EXE "${SCRIPT_DIR}/analysis_intervening_genes.py" -c "$CHROM_PICKLE" -sy "$SYNT_FOR_ANALYSIS" -ms "$CLUSTERS_FOR_ANALYSIS" -o "$INTERVENING_GENES_OUTPUT_PREFIX_NODE_CURRENT_RUN"
                if [ -f "$INTERVENING_GENES_TSV_NODE_CURRENT_RUN" ]; then
                    echo "Step 6 (Intervening Genes Analysis) for Node '$CURRENT_NODE_NAME' finished. Output: $INTERVENING_GENES_TSV_NODE_CURRENT_RUN";
                else
                    echo "WARNING: Step 6 output file ('$INTERVENING_GENES_TSV_NODE_CURRENT_RUN') not found for '$CURRENT_NODE_NAME'.";
                fi
            else
                echo "ERROR: Inputs for Step 6 (Intervening Genes Analysis) missing for '$CURRENT_NODE_NAME'."
                echo "       Needed .synt file: '$SYNT_FOR_ANALYSIS' (Exists? $(test -f $SYNT_FOR_ANALYSIS && echo Yes || echo No))"
                echo "       Needed .clusters file: '$CLUSTERS_FOR_ANALYSIS' (Exists? $(test -f $CLUSTERS_FOR_ANALYSIS && echo Yes || echo No))"
            fi
        else # RUN_INTERVENING_GENES_ANALYSIS is false
            echo "Skipping Step 6 (Optional Intervening Genes Analysis) for Node '$CURRENT_NODE_NAME' as RUN_INTERVENING_GENES_ANALYSIS is false."
        fi
    else # START_STEP > 6
        echo "Skipping Step 6 execution for Node '$CURRENT_NODE_NAME' (START_STEP=$START_STEP)."
    fi

    echo "Finished processing Node: $CURRENT_NODE_NAME"
    echo "###################################################################################"
    echo ""
done

SCRIPT_EXECUTION_END_DATETIME_ACTUAL=$(date -u +"%Y-%m-%d %H:%M:%S UTC")
echo "==================================================================================="
echo "SYNPHONI Multi-Node Pipeline finished for all specified nodes!"
echo "Execution End Time (Actual): $SCRIPT_EXECUTION_END_DATETIME_ACTUAL"
echo "Outputs are organized in subdirectories within: $CURRENT_RUN_BASE_OUTPUT_DIR"
echo "Please review log messages for warnings or errors for each node."
echo "==================================================================================="
