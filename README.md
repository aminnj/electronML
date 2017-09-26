# electronML

## Installing
```bash
CMSSW_VERSION=CMSSW_8_0_28
cmsrel $CMSSW_VERSION
cd $CMSSW_VERSION/src
git clone https://github.com/aminnj/electronML
cd electronML
scram b -j10
```

## Running
Use the pset in `python/run_cfg.py`
