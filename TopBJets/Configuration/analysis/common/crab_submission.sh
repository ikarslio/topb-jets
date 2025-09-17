rm -r data && rm -r data_v0 && rm -r data_v1 && rm -r data_v2 && rm -r data_v3 && rm -r data_v4
mkdir data && cd data && crab3_submit.py -p ../test/crab3/production/2017_106X_DoubleEG.json --Tier2 T2_US_Purdue
crab3_localarea.py --mask-data
mkdir data && cd data && crab3_submit.py -p ../test/crab3/production/2017_106X_MuonEG.json --Tier2 T2_US_Purdue
sed "s/'_OLD'/'_OLD_v0'/g" test/crab3/utilities/crab3_localarea.py test/crab3/utilities/crab3_localarea.py
crab3_localarea.py --mask-data
mkdir data && cd data && crab3_submit.py -p ../test/crab3/production/2017_106X_DoubleMuon.json --Tier2 T2_US_Purdue
sed "s/'_OLD_v0'/'_OLD_v1'/g" test/crab3/utilities/crab3_localarea.py test/crab3/utilities/crab3_localarea.py
crab3_localarea.py --mask-data
mkdir data && cd data && crab3_submit.py -p ../test/crab3/production/2017_106X_SingleElectron.json --Tier2 T2_US_Purdue
sed "s/'_OLD_v1'/'_OLD_v2'/g" test/crab3/utilities/crab3_localarea.py test/crab3/utilities/crab3_localarea.py
crab3_localarea.py --mask-data
mkdir data && cd data && crab3_submit.py -p ../test/crab3/production/2017_106X_SingleMuon.json --Tier2 T2_US_Purdue
sed "s/'_OLD_v2'/'_OLD_v3'/g" test/crab3/utilities/crab3_localarea.py test/crab3/utilities/crab3_localarea.py
crab3_localarea.py --mask-data
mkdir data && cd data && crab3_submit.py -p ../test/crab3/production/2018_106X_EGamma.json --Tier2 T2_US_Purdue
crab3_submit.py -p ../test/crab3/production/2018_106X_MuonEG.json --Tier2 T2_US_Purdue
sed "s/'_OLD_v3'/'_OLD_v4'/g" test/crab3/utilities/crab3_localarea.py test/crab3/utilities/crab3_localarea.py
mkdir data && cd data && crab3_submit.py -p ../test/crab3/production/2018_106X_DoubleMuon.json --Tier2 T2_US_Purdue
crab3_submit.py -p ../test/crab3/production/2018_106X_SingleMuon.json --Tier2 T2_US_Purdue
