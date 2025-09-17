Analysis scripts
================

* Apply pre-selection to 2017 and 2018 TopBJets ntuples
  * source run2017/analysisminitree.sh
  * source run2018/analysisminitree.sh

* Perform reconstruction on mini-Trees and measure signal yields
  * source analysiscripthist.sh

* Perform hadd
  * ```hadd -f data.root run2017/DoubleEG/ee_run2017.root run2017/MuonEG/emu_run2017.root run2017/DoubleMuon/mumu_run2017.root run2017/SingleElectron/se_run2017.root run2017/SingleMuon/smu_run2017.root run2018/EGamma/ee_run2018.root run2018/MuonEG/emu_run2018.root run2018/DoubleMuon/mumu_run2018.root run2018/SingleMuon/smu_run2018.root```

* Plotting
  * `root -l MassFitsData.C`
