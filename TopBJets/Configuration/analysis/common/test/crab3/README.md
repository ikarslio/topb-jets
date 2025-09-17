crab3 setup
===============

* `source crab3_env.sh`

* if starting a new production, `crab3_localarea.py --mask-data` (see script description below).

* prepare the **datasets JSON**, i.e. a JSON file containing the configuration of the TopBJets ntuples to be produced via crab3. This file needs to be prepared by hand. Examples can be found in the `datasets/` directory.

* run `crab3_production_json.py` to create the **production JSON**. This contains all the information needed for submitting a set of crab3 tasks. Necessary inputs include the datasets JSON(s), the target Tier-2 storage area and the target output directory for the final TopBJets ntuples. This JSON file should be kept at hand as long as the corresponding tasks are running (see `crab3_monitor.py`). For more info, do `crab3_production_json.py -h`.

  * `source json_production.sh`

* run `crab3_submit.py` to submit crab3 tasks. Each task is defined as an entry in the **production JSON** used as input to crab3_submit.py. For more info, do `crab3_submit.py -h`.

  * `crab3_submit.py -p production/2017_106X_DoubleEG.json --Tier2 T2_US_Purdue`

* run `crab3_monitor.py` periodically to check the status of the specified crab3 tasks. Option `--resubmit`: failed jobs will be resubmitted. Option `--hadd`: the outputs of a fully completed task are merged into the final TopBJets ntuple. Some options, e.g. --hadd, require the production JSON as input. For more info, do `crab3_monitor.py -h`.

  * `crab3_monitor.py -t crab_ee_run2017B --hadd -p production/2017_106X_DoubleEG.json --Tier2-prepath /eos/purdue --Tier2-prefix davs://eos.cms.rcac.purdue.edu:9000`

* when re-logging in to check on the tasks, do `source crab3_env.sh`, then run `crab3_monitor.py -t crab_ee_run2017B`.

Utilities
===============================

**utilities/** directory: contains scripts to facilitate simple routine tasks

* `crab3_env.sh`:

  * source this file to set up CMSSW and crab3 environments

* `crab3_localarea.py --mask-data`:

  * WARNING: executing this script will rename some directories, double-check before executing
  * renames `data/` sub-directories in `$CMSSW_BASE` that are not required for ntuple production
  * necessary to bring size of crab3 tarball below max-size limit
  * this script is supposed to be executed once, before submitting jobs
  * if not executed, grid jobs will most likely fail due to crab3 tarball exceeding max-size limit
