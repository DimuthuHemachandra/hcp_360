#!/bin/bash



in=$1 
out=$2 
analysis_level=$3
participant_label=$4

/neurodocker/startup.sh python workflow.py $in $out $analysis_level --participant_label $participant_label



