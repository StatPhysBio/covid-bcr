WDPATH=/gscratch/stf/zachmon/igor_wd/hiv/ #Let's define a shorthand for the working directory
MYCOMMANDS="igor -set_wd $WDPATH"
BATCHNAME=hiv_gen
MODELPATH=/gscratch/stf/zachmon/igor_wd/hiv_not_low_enough/hiv_inference/
PARMSFILE=${MODELPATH}final_parms.txt
MARGFILE=${MODELPATH}final_marginals.txt
MYCOMMANDS="$MYCOMMANDS -species human -chain heavy_naive -set_custom_model ${PARMSFILE} ${MARGFILE}" #Add chain and species commands
$MYCOMMANDS -batch $BATCHNAME -generate 1000000 #Generate
