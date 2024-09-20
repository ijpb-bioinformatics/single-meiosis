#! /bin/bash

# Date: 23/10/2017
# Author: Delphine Charif
# Script to process one Tetrad


# EX: ./singlemeiosis.sh tetrad_name V3 /work/dcharif/V3 delphine.charif@inra.fr /save/dcharif/bioinfo-dev/sh-dev/single-meiosis/share/singlemeiosis-pipeline/etc/singlemeiosis-pipeline_genotoul_V3

# Tetrade name
TETRAD=$1

# Analysis version
ANALYSIS_VERSION=$2

# Out data dir
OUT=$3

# mail
MAIL=$4

# single-meiosis pipeline Config file
CONFIG_FILE=$5

## Environment variables

# CODE

SH_CODE_sm=/sh-dev/single-meiosis/bin/singlemeiosis-pipeline.sh
CONFIG_FILE_sm=$CONFIG_FILE
R_CODE_Vutils=/r-dev/variantutils/doc/
R_CODE_HMM=/r-dev/hmm-nco/functions

# SCRIPTS soumission cluster

submitR=/sh-dev/common-scripts/cluster-submission-scripts/bin/submitR
submitSH=/bioinfo-code/bin/SGEqsub_v1 

# RAW DATA DIRECTORY

DATA_DIR=/raw-data/single-meiosis

# OUTPUTDIR

OUTDIR=$OUT/$ANALYSIS_VERSION
OUTDIR_TETRAD=$OUTDIR/$TETRAD
INPUTDIR_HMM=$OUTDIR_TETRAD/Analysis/Gold_Variant/InputHMM
OUTDIR_HMM=$OUTDIR_TETRAD/Analysis/Gold_Variant/OutputHMM

# Variant ressources

VDir=/save/project/ijpb/bioinfo-projects/single-meiosis/Col-Ler-Polymorphisms/SI_Schneeberger/Gold_Variant 
VPath=/save/project/ijpb/bioinfo-projects/single-meiosis/Col-Ler-Polymorphisms/SI_Schneeberger/Gold_Variant/newsdi_Gold_Variant.txt
VColPos=/save/project/ijpb/bioinfo-projects/single-meiosis/Col-Ler-Polymorphisms/SI_Schneeberger/Gold_Variant/dfColPos
VLerPos=/save/project/ijpb/bioinfo-projects/single-meiosis/Col-Ler-Polymorphisms/SI_Schneeberger/Gold_Variant/dfLerPos

## 1) Lancement du pipeline shell singlemeiosis-pipeline.sh

### 1-a) Generer le script de configuration du pipeline shell singlemeiosis-pipeline.sh pour la tétrade

cp $CONFIG_FILE_sm $OUTDIR/singlemeiosis-pipeline_$TETRAD.config

echo "
#[tetrad_samples]

# M1
m1_sample_name_alias=M1
m1_sample_seqfile_R1=$DATA_DIR/$TETRAD/$TETRAD\_M1_s1.fq.gz
m1_sample_seqfile_R2=$DATA_DIR/$TETRAD/$TETRAD\_M1_s2.fq.gz
# M2
m2_sample_name_alias=M2
m2_sample_seqfile_R1=$DATA_DIR/$TETRAD/$TETRAD\_M2_s1.fq.gz
m2_sample_seqfile_R2=$DATA_DIR/$TETRAD/$TETRAD\_M2_s2.fq.gz
# M3
m3_sample_name_alias=M3
m3_sample_seqfile_R1=$DATA_DIR/$TETRAD/$TETRAD\_M3_s1.fq.gz
m3_sample_seqfile_R2=$DATA_DIR/$TETRAD/$TETRAD\_M3_s2.fq.gz
# M4
m4_sample_name_alias=M4
m4_sample_seqfile_R1=$DATA_DIR/$TETRAD/$TETRAD\_M4_s1.fq.gz 
m4_sample_seqfile_R2=$DATA_DIR/$TETRAD/$TETRAD\_M4_s2.fq.gz
" >> $OUTDIR/singlemeiosis-pipeline_$TETRAD.config

### 1-b) Lancement du pipeline

jid=`bash $submitSH -s "$SH_CODE_sm -c $OUTDIR/singlemeiosis-pipeline_$TETRAD.config -o $OUTDIR/$TETRAD" -j $TETRAD-sm -M $MAIL -o $OUTDIR`
id=`echo $jid | awk '{print $3}'`
echo $jid

## 2) Lancement du script R Vutils pour la tétrade

jid1=`bash $submitR -j ${TETRAD}-Analysis -M $MAIL -w $id -c -n $OUTDIR $R_CODE_Vutils/Launch_Tetrad_Analysis_genotoul_V3.R $OUTDIR/$TETRAD $VDir $VPath $VColPos $VLerPos`
id1=`echo $jid1 | awk '{print $3}'`
echo $jid1

## 3) Lancement des scripts R HMM-nco

jid2=`bash $submitR -j ${TETRAD}-JointTetrade -M $MAIL -w $id1 -c -n $OUTDIR_HMM $R_CODE_HMM/FitContTimeMC-BB-JointTetrade.R $INPUTDIR_HMM $OUTDIR_HMM $TETRAD`
id2=`echo $jid2 | awk '{print $3}'`
echo $jid2

jid3=`bash  $submitR -j ${TETRAD}-MergeResults -M $MAIL -c -w $id2 -n $OUTDIR_HMM $R_CODE_HMM/MergeResultsHMM-BB-Tetrade.R $INPUTDIR_HMM $OUTDIR_HMM $TETRAD`
id3=`echo $jid3 | awk '{print $3}'`
echo $jid3

jid4=`bash $submitR -j ${TETRAD}-PlotResults -M $MAIL -c -w $id3 -n $OUTDIR_HMM $R_CODE_HMM/MergeResultsHMM-BB-Tetrade-Plot.R $INPUTDIR_HMM $OUTDIR_HMM $TETRAD`
echo $jid4


