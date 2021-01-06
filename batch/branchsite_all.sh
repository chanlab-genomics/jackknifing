#!/bin/bash
#PBS -l ncpus=1
#PBS -l mem=4GB
#PBS -j oe
#PBS -P d85

set -o errexit
cd $PBS_O_WORKDIR
source ~/miniconda3/bin/activate
conda activate ete_env

export ALG=/scratch/d85/kd0820/selection/codon
export TRE=/scratch/d85/kd0820/selection/trees_branch_site
export CTL=/scratch/d85/kd0820/selection/control_files
export IDS=/scratch/d85/kd0820/selection/orthogroup_lists

export GRP=/scratch/d85/kd0820/selection/BranchSite_Raw_out/$ORTHOLIST/${ORTHOLIST}_${MODEL}
export OUT=/scratch/d85/kd0820/selection/BranchSite_Sub_out/$ORTHOLIST/${ORTHOLIST}_${MODEL}/out
export LNL=/scratch/d85/kd0820/selection/BranchSite_Sub_out/$ORTHOLIST/${ORTHOLIST}_${MODEL}/lnl

#############
#############

# Go to directory containing all raw output from codeml analysis
cd $GRP
ete3 evol -v 3 -C 1 --noimg --codeml_binary ~/apps/paml-v4.9j/bin/codeml --slr_binary ~/apps/Slr \
        -t $TRE/$ORTHONAME.branch_site.nw \
        --codeml_config_file $CTL/${MODEL}_codeml.ctl \
        --alg $ALG/$ORTHONAME.edited.codon \
        -o ${JOBNAME}

# copy out file to separate directory for text parsing results
cp $JOBNAME/XX*/out $OUT/${JOBNAME}_out
sed -n '/^lnL/p' $OUT/${JOBNAME}_out > $LNL/${JOBNAME}_out_lnL
sed -n '/^kappa /p' $OUT/${JOBNAME}_out > $LNL/${JOBNAME}_out_kappa
sed -n '/^proportion/p' $OUT/${JOBNAME}_out > $LNL/${JOBNAME}_out_proportion
sed -n '/^background w/p' $OUT/${JOBNAME}_out > $LNL/${JOBNAME}_out_background
sed -n '/^foreground w/p' $OUT/${JOBNAME}_out > $LNL/${JOBNAME}_out_foreground

# compress raw output directory
tar -cvf $JOBNAME.tar $JOBNAME/
if [ $? != 0 ]
then
 echo "Backup FAILED for ${JOBNAME}"
else
 echo "Backup SUCCESS for ${JOBNAME}"
 rm -rf $JOBNAME
fi


# move stdout files for successful codeml runs
cd $PBS_O_WORKDIR
bash check_successful_completion_BS.sh
