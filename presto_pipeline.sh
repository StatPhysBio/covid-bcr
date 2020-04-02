#!/bin/bash
#SBATCH --job-name=presto
#SBATCH -p stf
#SBATCH -A stf
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=28
#SBATCH --time=5-00:00:00

#INPUTDIR=${2}
SCRATCHDIR=/gscratch/stf/zachmon/
SOFTWAREDIR=${SCRATCHDIR}software/presto-0.5.13/bin/
COVIDDIR=${SCRATCHDIR}covid/
ALIGNDIR=${COVIDDIR}aligned/
FILTERDIR=${COVIDDIR}filtered/
TRIMDIR=${COVIDDIR}masked/
DEDUPDIR=${COVIDDIR}collapsed/
PRIMERDIR=${COVIDDIR}data/primers/

usage() { echo "Usage: $0 [-s 1-18][-a] [-f] [-m] [-c]" 1>&2; exit 1; }
while getopts ":s:afmc" arg; do
    case $arg in
        s)
            SAMPLENUM=${OPTARG}
            INPUTDIR=${COVIDDIR}data/sample_${SAMPLENUM}/
            FILENAME=S${SAMPLENUM}
            ;;
        a)
            ASSEMBLE=true
            ;;
        f)
            FILTER=true
            ;;
        m)
            MASK=true
            ;;
        c)
            COLLAPSE=true
            ;;
        *) 
            usage
            ;;
    esac
done
shift $((OPTIND-1))

echo "SAMPLE NUM $SAMPLENUM"
echo "ASSEMBLE $ASSEMBLE"
echo "FILTER $FILTER"
echo "MASK $MASK"
echo "COLLAPSE $COLLAPSE"
#if [ -z "${s}" ] || [ -z "${a}" ] || [ -z "${f}" ] || [ -z "${m}" ] || [ -z "${c}" ]; then
#    usage
#fi
#========== Assembling Pairs ==========
if [ "$ASSEMBLE" = true ]; then
    if [ ! -d "$ALIGNDIR" ]; then
      mkdir $ALIGNDIR
    fi

    #  Assemble pairs.
    ${SOFTWAREDIR}AssemblePairs.py align -1 ${INPUTDIR}${SAMPLENUM}-IGG_1.fastq \
                                         -2 ${INPUTDIR}${SAMPLENUM}-IGG_2.fastq \
                                         --coord illumina \
                                         --rc tail \
                                         --outname ${ALIGNDIR}${FILENAME} \
                                         --log ${ALIGNDIR}assemble_S${SAMPLENUM}.log
fi
#========== Filter Sequences ==========
if [ "$FILTER" = true ]; then
    if [ ! -d "$FILTERDIR" ]; then
      mkdir $FILTERDIR
    fi

    #  Filter sequences using a QScore of 30.
    ${SOFTWAREDIR}FilterSeq.py quality -s ${ALIGNDIR}${FILENAME}_assemble-pass.fastq \
                                       -q 30 \
                                       --outname ${FILTERDIR}30_${FILENAME} \
                                       --log ${FILTERDIR}filter_30_S${SAMPLENUM}.log
fi
#========== Trim/Cut Primers  ==========
if [ "$MASK" = true ]; then
    if [ ! -d "$TRIMDIR" ]; then
      mkdir $TRIMDIR
    fi

    #  Mask V primers.
    ${SOFTWAREDIR}MaskPrimers.py score -s ${FILTERDIR}30_${FILENAME}_quality-pass.fastq \
                                       -p ${PRIMERDIR}vprimers.fa \
                                       --mode mask --pf VPRIMER \
                                       --outname ${TRIMDIR}${FILENAME}-FWD \
                                       --log ${TRIMDIR}${SAMPLENUM}PV.log

    #  Cut C primer.
    ${SOFTWAREDIR}MaskPrimers.py score -s ${TRIMDIR}${FILENAME}-FWD_primers-pass.fastq \
                                       -p ${PRIMERDIR}cprimers.fa \
                                       --mode cut --revpr --pf CPRIMER \
                                       --outname ${FILENAME}-REV \
                                       --log ${TRIMDIR}${SAMPLENUM}PC.log
fi
#========== Deduplicate Sequences ==========
if [ "$COLLAPSE" = true ]; then
    if [ ! -d "$DEDUPDIR" ]; then
      mkdir $DEDUPDIR
    fi

    #  Deduplicate data.
    ${SOFTWAREDIR}CollapseSeq.py -s ${TRIMDIR}${FILENAME}-REV_primers-pass.fastq \
                                 -n 0 \
                                 --uf CPRIMER --cf VPRIMER --act set --inner \
                                 --outname ${DEDUPDIR}unique_${FILENAME}
fi
