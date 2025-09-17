find . -name \run2018*_data*.root -type f -delete
sed -i 's+20230523+20230731+g' */prepareslurmscript.sh
cd EGamma && source prepareslurmscript.sh && source RuncmsSlurm_slurmscriptrun.sh && cd ..
cp EGamma/outfiles1/MiniTree1.C MuonEG/outfiles1 && cd MuonEG && source prepareslurmscript.sh && source RuncmsSlurm_slurmscriptrun.sh && cd ..
cp EGamma/outfiles1/MiniTree1.C DoubleMuon/outfiles1 && cd DoubleMuon && source prepareslurmscript.sh && source RuncmsSlurm_slurmscriptrun.sh && cd ..
cp EGamma/outfiles1/MiniTree1.C SingleMuon/outfiles1 && cd SingleMuon && source prepareslurmscript.sh && source RuncmsSlurm_slurmscriptrun.sh && cd ..
