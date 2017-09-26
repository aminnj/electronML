# electronML

This is a modified version of the official ntupler used for the BDT-based classification of electrons
used [here](https://github.com/Werbellin/ID_flat_ntupler/blob/master/Ntuplizer/python/run_AOD_cfg.py).

## Installing
```bash
CMSSW_VERSION=CMSSW_8_0_28
cmsrel $CMSSW_VERSION
cd $CMSSW_VERSION/src
mkdir Ntuplizer
cd Ntuplizer
git clone https://github.com/aminnj/electronML
cd ..
scram b -j10
cd Ntuplizer/electronML
```

## Running
Use the pset in `python/run_cfg.py`
