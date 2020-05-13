#!/bin/bash
#SBATCH --job-name=covid_preprocess
#SBATCH -p stf
#SBATCH -A stf
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=28
#SBATCH --time=5-00:00:00

# TODO Take in a working directory
#      Take in a conda directory for error correction
SCRATCHDIR=/gscratch/stf/zachmon/
CONDADIR=${SCRATCHDIR}miniconda3/
CONDAENVDIR=${CONDADIR}envs/covid-bcr
RUNCONDA=${CONDADIR}etc/profile.d/conda.sh
SOFTWAREDIR=${SCRATCHDIR}software/presto-0.5.13/bin/
COVIDDIR=${SCRATCHDIR}covid/round2data/
ALIGNDIR=${COVIDDIR}aligned/
FILTERDIR=${COVIDDIR}filtered/
TRIMDIR=${COVIDDIR}masked/
DEDUPDIR=${COVIDDIR}collapsed/
PRIMERDIR=${COVIDDIR}raw_data/primers/
DEDUPONELINE=${COVIDDIR}oneline_collapsed/
ERRORDIR=${COVIDDIR}error_corrected/

usage() { echo "Usage: $0 [-s <string>] [-a] [-f] [-m] [-c] ][-e] [-h]" 1>&2; exit 1; }
while getopts ":s:afmceh" arg; do
    case $arg in
        s) SAMPLENUM=${OPTARG}
           INPUTDIR=${COVIDDIR}raw_data/sample_${SAMPLENUM}/
           FILENAME=S${SAMPLENUM}
           echo "Input dir is ${INPUTDIR}";;
        a) ASSEMBLE=true;;
        f) FILTER=true;;
        m) MASK=true;;
        c) COLLAPSE=true;;
        e) ERRORCORRECT=true;;
        h)  usage;  exit 0;;
        \? ) echo "Unknown option: -$OPTARG\n" >&2; exit 1;;
        : ) echo -e "Missing argument for -$OPTARG\n" >&2; exit 1;;
        * ) echo -e "Unimplemented option: -$OPTARG\n" >&2; exit 1;;
    esac
done
shift $((OPTIND-1))

#  Ensure sample num is an integer.
#if ! [[ "$SAMPLENUM" =~ ^[0-9]+$ ]]; then
#    echo "Sample num must be an interger."
#    usage
#    exit 1
#fi

#  Ensure sample num is between 1 and 18.
#  TODO Update when more samples become available.
#if [[ "$SAMPLENUM" -gt 18 ]]; then
#    echo "Sample num is between 1 and 18." >&2
#    usage
#    exit 1
#fi
echo "SAMPLE NUM $SAMPLENUM"
echo "ASSEMBLE ${ASSEMBLE}"
echo "FILTER $FILTER"
echo "MASK $MASK"
echo "COLLAPSE $COLLAPSE"
echo "ERROR CORRECT $ERRORCORRECT"

#if [ -z "${a}" ] && [ -z "${f}" ] && [ -z "${m}" ] && [ -z "${c}" ]; then
#    usage
#fi

#========== Assembling Pairs ==========
if [ "${ASSEMBLE}" = true ]; then
    if [ ! -d "$ALIGNDIR" ]; then
      mkdir $ALIGNDIR
    fi
    echo "ASSEMBLING!"
    #  Assemble pairs.
    ${SOFTWAREDIR}AssemblePairs.py align -1 ${INPUTDIR}IgG${SAMPLENUM}_1.fastq \
                                         -2 ${INPUTDIR}IgG${SAMPLENUM}_2.fastq \
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
    echo "FILTERING!"
    #  Filter sequences using a QScore of 30.
    ${SOFTWAREDIR}FilterSeq.py quality -s ${ALIGNDIR}${FILENAME}_assemble-pass.fastq \
                                       --fasta \
                                       -q 30 \
                                       --outname ${FILTERDIR}30_${FILENAME} \
                                       --log ${FILTERDIR}filter_30_S${SAMPLENUM}.log
fi
#========== Trim/Cut Primers  ==========
if [ "$MASK" = true ]; then
    if [ ! -d "$TRIMDIR" ]; then
      mkdir $TRIMDIR
    fi
    echo "MASKING!"
    #  Mask V primers.
    ${SOFTWAREDIR}MaskPrimers.py score -s ${FILTERDIR}30_${FILENAME}_quality-pass.fasta \
                                       -p ${PRIMERDIR}vprimers.fa \
                                       --mode mask --pf VPRIMER \
                                       --outname ${TRIMDIR}${FILENAME}-FWD \
                                       --log ${TRIMDIR}${SAMPLENUM}PV.log

    #  Cut C primer.
    ${SOFTWAREDIR}MaskPrimers.py score -s ${TRIMDIR}${FILENAME}-FWD_primers-pass.fasta \
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
    if [ ! -d "$DEDUPONELINE" ]; then
      mkdir $DEDUPONELINE
    fi
    echo "DEDUPLICATING!"
    #  Deduplicate data.
    ${SOFTWAREDIR}CollapseSeq.py -s ${TRIMDIR}${FILENAME}-REV_primers-pass.fasta \
                                 -n 20 \
                                 --uf CPRIMER --cf VPRIMER --act set --inner \
                                 --outname ${DEDUPDIR}unique_${FILENAME}
   
    #  Convert sequence in fasta file to one line.
    fasta_formatter -i ${DEDUPDIR}unique_${FILENAME}_collapse-unique.fasta -o ${DEDUPONELINE}${FILENAME}_collapse-unique.fasta
fi

#========== Error Correct Singleton Sequences ==========
if [ "$ERRORCORRECT" = true ]; then
    if [ ! -d "$ERRORDIR" ]; then
      mkdir $ERRORDIR
    fi

    #  Export conda functions to be available to subshell.
    source ${RUNCONDA} 
    
    #  Create conda environment if necessary.
    if [ ! -d "${CONDAENVDIR}" ]; then
      echo "Conda environment covid-bcr does not exist. Creating now."
      conda env create -f env.yml
    fi
    
    #  Change to conda environment.
    conda activate covid-bcr
    
    #  Error correct data.
    python error_correct.py ${DEDUPONELINE}${FILENAME}_collapse-unique.fasta > ${ERRORDIR}${FILENAME}_corrected.fasta
    
    #  Deactivate conda environment.
    conda deactivate
fi
