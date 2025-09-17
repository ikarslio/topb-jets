#!/bin/sh
source /cvmfs/oasis.opensciencegrid.org/osg-software/osg-wn-client/3.6/current/el7-x86_64/setup.sh
export X509_USER_PROXY=~/x509up_u`id -u`
voms-proxy-init -voms cms -rfc -valid 192:00:00
sbatch cmsSlurmJobs/SlurmJob_1.sh
sbatch cmsSlurmJobs/SlurmJob_2.sh
sbatch cmsSlurmJobs/SlurmJob_3.sh
sbatch cmsSlurmJobs/SlurmJob_4.sh
sbatch cmsSlurmJobs/SlurmJob_5.sh
sbatch cmsSlurmJobs/SlurmJob_6.sh
sbatch cmsSlurmJobs/SlurmJob_7.sh
sbatch cmsSlurmJobs/SlurmJob_8.sh
sbatch cmsSlurmJobs/SlurmJob_9.sh
sbatch cmsSlurmJobs/SlurmJob_10.sh
sbatch cmsSlurmJobs/SlurmJob_11.sh
sbatch cmsSlurmJobs/SlurmJob_12.sh
sbatch cmsSlurmJobs/SlurmJob_13.sh
sbatch cmsSlurmJobs/SlurmJob_14.sh
sbatch cmsSlurmJobs/SlurmJob_15.sh
sbatch cmsSlurmJobs/SlurmJob_16.sh
sbatch cmsSlurmJobs/SlurmJob_17.sh
sbatch cmsSlurmJobs/SlurmJob_18.sh
sbatch cmsSlurmJobs/SlurmJob_19.sh
sbatch cmsSlurmJobs/SlurmJob_20.sh
sbatch cmsSlurmJobs/SlurmJob_21.sh
sbatch cmsSlurmJobs/SlurmJob_22.sh
sbatch cmsSlurmJobs/SlurmJob_23.sh
sbatch cmsSlurmJobs/SlurmJob_24.sh
