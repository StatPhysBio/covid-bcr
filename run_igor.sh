#!/bin/bash
INPUTDIR=/gscratch/stf/zachmon/covid/igor_input/
mapfile -t fasta_files < <( ls $INPUTDIR | awk '{ print $1 }' )
for j in ${fasta_files[@]}
do
    SUBSTRING=${j%.*}
    SUFFIX=${j##*.}
    COHORT=${j%%_*}    
    if [[ "$COHORT" != "all" ]]; then
        continue
    fi
    echo $SUBSTRING
    WD=/gscratch/stf/zachmon/covid/igor_wd/${SUBSTRING}_gen
    FASTA=${INPUTDIR}${j}
    sbatch igor_script.sh ${WD} ${FASTA} ${COHORT} 
done
