#!/bin/bash
#SBATCH --job-name=igor
#SBATCH -p spe 
#SBATCH -A spe
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=40
#SBATCH --time=6:00:00
#IGOR13="/gscratch/stf/zachmon/software/igor_1-3-0_exec/igor"
IGOR14="/gscratch/stf/zachmon/software/igor_1-4-0_exec/igor"
WDPATH=${1}
FASTA=${2}
BATCHNAME=${3}
ABSTARV="/gscratch/stf/zachmon/covid/covid-bcr/two_allele_abstar_ungapped_IGH.fasta"
ANCHORV="/gscratch/stf/zachmon/covid/covid-bcr/igor_V_gene_CDR3_anchors.csv"
ABSTARJ="/gscratch/stf/zachmon/covid/covid-bcr/abstar_genomic_Js.fasta"
ABSTARD="/gscratch/stf/zachmon/covid/covid-bcr/abstar_genomic_Ds.fasta"
ANCHORJ="/gscratch/stf/zachmon/covid/covid-bcr/igor_J_gene_CDR3_anchors.csv"
echo "WDPATH ${WDPATH}"
echo "FASTA ${FASTA}"
echo "BATCHNAME ${BATCHNAME}"
rm -rf $WDPATH
mkdir $WDPATH
MYCOMMANDS="${IGOR14} -set_wd ${WDPATH} -threads 40"
${MYCOMMANDS} -batch ${BATCHNAME} -read_seqs ${FASTA} #Read seqs
MYCOMMANDS="$MYCOMMANDS -species human -chain heavy_naive" #Add chain and species commands
MYCOMMANDS="$MYCOMMANDS -set_genomic --V ${ABSTARV} -set_CDR3_anchors --V ${ANCHORV}" #  Set v genes and anchors for abstar
MYCOMMANDS="$MYCOMMANDS -set_genomic --J ${ABSTARJ} -set_CDR3_anchors --J ${ANCHORJ}" #  Set j genes and anchors for abstar
MYCOMMANDS="$MYCOMMANDS -set_genomic --D ${ABSTARD}" #  Set d genes
${MYCOMMANDS} -batch ${BATCHNAME} -align --all #Align
${MYCOMMANDS} -batch ${BATCHNAME} -infer --L_thresh "1e-300" --N_iter 10 #  Infer
${MYCOMMANDS} -batch ${BATCHNAME} -evaluate -output --scenarios 10 #Evaluate
${MYCOMMANDS} -batch ${BATCHNAME} -generate 10 #Generate
