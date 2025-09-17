# Remove 'rm slurm*.out' from the script when running for the first time
cd run2017 &&
cd DoubleEG && find . -name \*data*.root > outputFileList.txt && rm slurm*.out && sbatch slurmscripthist.sh && cd .. &&
cp DoubleEG/DataYield.C DoubleEG/slurmscripthist.sh MuonEG && cd MuonEG && sed -i 's+ee_run2017+emu_run2017+g' DataYield.C && find . -name \*data*.root > outputFileList.txt && rm slurm*.out && sbatch slurmscripthist.sh && cd .. &&
cp DoubleEG/DataYield.C DoubleEG/slurmscripthist.sh DoubleMuon && cd DoubleMuon && sed -i 's+ee_run2017+mumu_run2017+g' DataYield.C && find . -name \*data*.root > outputFileList.txt && rm slurm*.out && sbatch slurmscripthist.sh && cd .. &&
cp DoubleEG/DataYield.C DoubleEG/slurmscripthist.sh SingleElectron && cd SingleElectron && sed -i 's+ee_run2017+se_run2017+g' DataYield.C && find . -name \*data*.root > outputFileList.txt && rm slurm*.out && sbatch slurmscripthist.sh && cd .. &&
cp DoubleEG/DataYield.C DoubleEG/slurmscripthist.sh SingleMuon && cd SingleMuon && sed -i 's+ee_run2017+smu_run2017+g' DataYield.C && find . -name \*data*.root > outputFileList.txt && rm slurm*.out && sbatch slurmscripthist.sh && cd ../.. &&
cp run2017/DoubleEG/DataYield.C run2017/DoubleEG/slurmscripthist.sh run2018/EGamma &&
cd run2018 &&
cd EGamma && sed -i 's+ee_run2017+ee_run2018+g' DataYield.C && find . -name \*data*.root > outputFileList.txt && rm slurm*.out && sbatch slurmscripthist.sh && cd .. &&
cp EGamma/DataYield.C EGamma/slurmscripthist.sh MuonEG && cd MuonEG && sed -i 's+ee_run2018+emu_run2018+g' DataYield.C && find . -name \*data*.root > outputFileList.txt && rm slurm*.out && sbatch slurmscripthist.sh && cd .. &&
cp EGamma/DataYield.C EGamma/slurmscripthist.sh DoubleMuon && cd DoubleMuon && sed -i 's+ee_run2018+mumu_run2018+g' DataYield.C && find . -name \*data*.root > outputFileList.txt && rm slurm*.out && sbatch slurmscripthist.sh && cd .. &&
cp EGamma/DataYield.C EGamma/slurmscripthist.sh SingleMuon && cd SingleMuon && sed -i 's+ee_run2018+smu_run2018+g' DataYield.C && find . -name \*data*.root > outputFileList.txt && rm slurm*.out && sbatch slurmscripthist.sh && cd ../..
#hadd -f data.root run2017/DoubleEG/ee_run2017.root run2017/MuonEG/emu_run2017.root run2017/DoubleMuon/mumu_run2017.root run2017/SingleElectron/se_run2017.root run2017/SingleMuon/smu_run2017.root run2018/EGamma/ee_run2018.root run2018/MuonEG/emu_run2018.root run2018/DoubleMuon/mumu_run2018.root run2018/SingleMuon/smu_run2018.root
