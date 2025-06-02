#!/bin/bash

# List of required directories
required_dirs=("../fastaq" "./ref")

# Check existence of each directory
missing_dirs=()
for dir in "${required_dirs[@]}"; do
    if [ ! -d "$dir" ]; then
        missing_dirs+=("$dir")
    fi
done

# If there are missing directories, display error message and exit
if [ ${#missing_dirs[@]} -ne 0 ]; then
    echo "Error: The following required directories are missing:"
    for dir in "${missing_dirs[@]}"; do
        echo "  - $dir"
    done
    echo "Process aborted."
    exit 1
fi

# Continue processing if all directories exist
echo "Required directory check completed. Proceeding with process."
echo "----------------------------------------"

# List of scripts
scripts=(
  "0010ref_files_processor.py"
  "0020create_list.py"
  "0110bbduk.py"
  "0120stroboalign.py"
  "0130picard_markduplicates.py"
  "0140create_metrics.py"
  "0210freebayes.py"
  "0310Beagle.py"
  "0320bcftools_consensus.py"
  "0330extract_vgsc_indels.py"
  "0410create_alignment.py"
  "0420extract_exons.py"
  "0430create_transcripts.py"
  "0510calculate_hamming_distance.py"
  "0520Create_Histogram_from_Distance_Matrix.py"
  "0610CDS_variant_analyzer.py"
  "0620remove_duplications.py"
  "0630add_genomeinfo.py"
  "0640merge_non_common_transcripts.py"
  "0650CDS_position_mapper.py"
  "0710create_DNAtable.py"
  "0720create_AAtable.py"
  "0730add_exon_name.py"
)

total_scripts=${#scripts[@]}
current_script=1

# Execute each script
for script in "${scripts[@]}"; do
    echo "Executing ($current_script/$total_scripts): $script"
    echo "----------------------------------------"

    if ! $script; then
        echo "Error: An error occurred while executing $script"
        exit 1
    fi

    echo "Completed: $script"
    echo "----------------------------------------"
    ((current_script++))
done

echo "All processes completed"