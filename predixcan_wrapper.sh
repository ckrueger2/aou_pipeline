#!/bin/bash

#command
usage() {
    echo "Usage: $0 --phecode <PHECODE> --pop <POP> --ref <REF> [--gwas_h2 <H2>] [--gwas_N <N>] [--databases]"
    exit 1
}
        
TISSUES="Adipose_Subcutaneous\nAdipose_Visceral_Omentum\nAdrenal_Gland\nArtery_Aorta\nArtery_Coronary\nArtery_Tibial\nBrain_Amygdala\nBrain_Anterior_cingulate_cortex_BA24\nBrain_Caudate_basal_ganglia\nBrain_Cerebellar_Hemisphere\nBrain_Cerebellum\nBrain_Cortex\nBrain_Frontal_Cortex_BA9\nBrain_Hippocampus\nBrain_Hypothalamus\nBrain_Nucleus_accumbens_basal_ganglia\nBrain_Putamen_basal_ganglia\nBrain_Spinal_cord_cervical_c-1\nBrain_Substantia_nigra\nBreast_Mammary_Tissue\nCells_Cultured_fibroblasts\nCells_EBV-transformed_lymphocytes\nColon_Sigmoid\nColon_Transverse\nEsophagus_Gastroesophageal_Junction\nEsophagus_Mucosa\nEsophagus_Muscularis\nHeart_Atrial_Appendage\nHeart_Left_Ventricle\nKidney_Cortex\nLiver\nLung\nMinor_Salivary_Gland\nMuscle_Skeletal\nNerve_Tibial\nOvary\nPancreas\nPituitary\nProstate\nSkin_Not_Sun_Exposed_Suprapubic\nSkin_Sun_Exposed_Lower_leg\nSmall_Intestine_Terminal_Ileum\nSpleen\nStomach\nTestis\nThyroid\nUterus\nVagina\nWhole_Blood"

#command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --phecode)
            PHECODE=$2
            shift 2
            ;;
        --pop)
            POP=$2
            shift 2
            ;;
        --ref)
            REF=$2
            shift 2
            ;;
        --gwas_h2)
            H2=$2
            shift 2
            ;;
        --gwas_N)
            N=$2
            shift 2
            ;;
        --databases)
            echo -e "$TISSUES"
            shift 1
            ;;
        *)
            echo "unknown flag: $1"
            usage
            ;;
    esac
done

#check for required arguments
if [[ -z "$PHECODE" || -z "$POP" || -z "$REF" ]]; then
    usage
fi

#github repo path
REPO=$HOME/GWAS-TWAS-in-All-of-Us-Cloud

#set up S-PrediXcan environment
bash "$REPO/set-up-predixcan.sh"

#activate conda
source ~/miniconda3/bin/activate

#create environment with compatible versions (version numbers may need to be changed with future updates)
if ! conda env list | grep -q imlabtools; then
    # If it doesn't exist, create it
    conda create -n imlabtools python=3.8 numpy pandas scipy -y
fi

#activate imlabtools
if conda activate imlabtools; then
    echo "Successfully activated imlabtools environment"
fi

#patch MetaXcan code
if [ -f /home/jupyter/MetaXcan/software/metax/gwas/GWAS.py ]; then
    sed -i 's/if a.dtype == numpy.object:/if a.dtype == object or str(a.dtype).startswith("object"):/' /home/jupyter/MetaXcan/software/metax/gwas/GWAS.py
fi

#patch numpy.str deprecation in Utilities.py
if [ -f /home/jupyter/MetaXcan/software/metax/metaxcan/Utilities.py ]; then
    # First check if the file contains the original numpy.str (not already patched)
    if grep -q "numpy\.str[^_]" /home/jupyter/MetaXcan/software/metax/metaxcan/Utilities.py; then
        # Only replace numpy.str with numpy.str_ if it hasn't been replaced yet
        sed -i 's/numpy\.str\([^_]\)/numpy.str_\1/g' /home/jupyter/MetaXcan/software/metax/metaxcan/Utilities.py
    elif grep -q "numpy\.str_" /home/jupyter/MetaXcan/software/metax/metaxcan/Utilities.py; then
        # If we find numpy.str__ (double underscore), replace it with numpy.str_
        sed -i 's/numpy\.str__/numpy.str_/g' /home/jupyter/MetaXcan/software/metax/metaxcan/Utilities.py
    fi
    
    # Fix the specific line that causes errors with numpy.str__
    sed -i 's/type = \[numpy\.str[_]*, numpy\.float64, numpy\.float64, numpy\.float64\]/type = \[str, numpy.float64, numpy.float64, numpy.float64\]/g' /home/jupyter/MetaXcan/software/metax/metaxcan/Utilities.py
    
    # Fix pandas drop() method call
    sed -i 's/results = results.drop("n_snps_in_model",1)/results = results.drop(columns=["n_snps_in_model"])/' /home/jupyter/MetaXcan/software/metax/metaxcan/Utilities.py
fi

output_file="/home/jupyter/${POP}_predixcan_output_${PHECODE}_${REF}.csv"

#check if the output file already exists
if [ -f "$output_file" ]; then
    echo "WARNING: Output file $output_file already exists."
    read -p "Press ENTER to replace it, or type 'n' to cancel: " response
    
    if [[ $response =~ ^[Nn]$ ]]; then
        echo "Operation cancelled by user."
        exit 1
    else
        # Delete the file
        rm -f "$output_file"
        echo "Existing file has been deleted."
    fi
fi

#run s-predixcan
PREDIXCAN_CMD="python $REPO/05run-predixcan.py --phecode \"$PHECODE\" --pop \"$POP\" --ref \"$REF\""
if [[ ! -z "$H2" ]]; then
    PREDIXCAN_CMD="$PREDIXCAN_CMD --gwas_h2 \"$H2\""
fi
if [[ ! -z "$N" ]]; then
    PREDIXCAN_CMD="$PREDIXCAN_CMD --gwas_N \"$N\""
fi
eval $PREDIXCAN_CMD

#run qqman on twas sum stats
Rscript "$REPO/06twas_qqman.R" --phecode "$PHECODE" --pop "$POP" --ref "$REF"

#deactivate imlabtools
conda deactivate

#how to view generated PNG files
echo "To view the S-PrediXcan and PNG files, go to the Jupyter file browser by selecting the jupyter logo to the top left of the terminal."