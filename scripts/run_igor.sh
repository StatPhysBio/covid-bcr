#!/bin/bash
#
#    Copyright (C) Montague, Zachary
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#    Script used to run IGoR to obtain the recombination model.
#    The user should change the variable so that they point to the correct
#    locations in their own directories.
#
#    Required software:
#    IGoR https://qmarcou.github.io/IGoR/ 

#SBATCH --job-name=igor
#SBATCH -p spe 
#SBATCH -A spe
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --ntasks-per-node=40
#SBATCH --time=6:00:00

IGOR14="/gscratch/stf/zachmon/software/igor_1-4-0_exec/igor"
WDPATH=/gscratch/stf/zachmon/covid/total_bcell/igor_wd/19_08_IGoR/
FASTA=/gscratch/stf/zachmon/covid/total_bcell/igor_input/abstar_all_individuals_most_common.fasta
BATCHNAME=foo
ABSTARV="/gscratch/stf/zachmon/covid/covid-bcr/igor_input/abstar_genomic_Vs_for_igor.fasta"
ANCHORV="/gscratch/stf/zachmon/covid/covid-bcr/igor_input/igor_V_gene_CDR3_anchors.csv"
ABSTARJ="/gscratch/stf/zachmon/covid/covid-bcr/igor_input/abstar_genomic_Js_for_igor.fasta"
ABSTARD="/gscratch/stf/zachmon/covid/covid-bcr/igor_input/abstar_genomic_Ds.fasta"
ANCHORJ="/gscratch/stf/zachmon/covid/covid-bcr/igor_input/igor_J_gene_CDR3_anchors.csv"
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
${MYCOMMANDS} -batch ${BATCHNAME} -align --all #  Align
${MYCOMMANDS} -batch ${BATCHNAME} -infer --L_thresh "1e-300" --N_iter 10 #  Infer
${MYCOMMANDS} -batch bar -generate 1000000 #  Generate
