#!/bin/bash
#SBATCH --job-name=igor
#SBATCH -p stf
#SBATCH -A stf
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=28
#SBATCH --time=2:00:00
#IGOR13="/gscratch/stf/zachmon/software/igor_1-3-0_exec/igor"
IGOR14="/gscratch/stf/zachmon/software/igor_1-4-0_exec/igor"
WDPATH=${1}
FASTA=${2}
BATCHNAME=${3}
echo "WDPATH ${WDPATH}"
echo "FASTA ${FASTA}"
echo "BATCHNAME ${BATCHNAME}"
rm -rf $WDPATH
mkdir $WDPATH
MYCOMMANDS="${IGOR14} -set_wd ${WDPATH} -threads 28"
${MYCOMMANDS} -batch ${BATCHNAME} -read_seqs ${FASTA} #Read seqs
MYCOMMANDS="$MYCOMMANDS -species human -chain heavy_naive" #Add chain and species commands
${MYCOMMANDS} -batch ${BATCHNAME} -align --all #Align
${MYCOMMANDS} -batch ${BATCHNAME} -infer --L_thresh "1e-300" --N_iter 10 #  Infer
${MYCOMMANDS} -batch ${BATCHNAME} -evaluate -output --scenarios 10 #Evaluate
${MYCOMMANDS} -batch ${BATCHNAME} -generate 10 #Generate
