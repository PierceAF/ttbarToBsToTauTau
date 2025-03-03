## Installation
- This analysis framework requires only ROOT and GCC (you don't need to set up CMSSW but need to check ROOT version dependency). You need singularity if you're system is OS is Alma9.
- For HiPerGator:
- ```source /cvmfs/cms.cern.ch/cmsset_default.sh```
- ```cmssw-el7```
- ```cmsrel CMSSW_10_2_16_UL```
- ```cd CMSSW_10_2_16_UL/src```
- ```cmsenv```
- ```git clone git@github.com:cghuh/ttbarToBsToTauTau.git```

Run outside of singularity
- ```python3 condor/setup.py```

After running the setup.py, the filelists are automatically generated using condor/filelist_2018.txt

### Run code locally
- Do it inside of cmssw-el7 singularity
- ```source setup.sh```
- ```make -j6 Analyzer Plotter```
- ```./Analyzer output.root filelists/2018/backgrounds/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.txt year=2018```
- additional option: ```ilast=NEVENT```, ```debug=True``` (default is False)


### Run using condor (Please run outside of singularity)
- ```python3 condor/run_all.py --full --nohadd --batch --condor --run --outdir=results/run_2024_11_28_syst --optim --nevt=1000000```
- Recover jobs when run_all.py is dead. ```python3 condor/run_all.py --full --nohadd --batch --condor --run --outdir=results/run_2024_11_28_syst --recover```


### Plotting after full condor run
- Do it inside of cmssw-el7 singularity
- When the condor jobs are done, you need to merge by hand(manually)
- ```./Plotter plot.root output.root```

## Explanation of code
- ```Analyzer.cc```: main code
- ```settings.h```: Set option (systematic on/off, scale factor on/off)
- ```include/Variables.h```: Define object and event variables
- ```include/Weighting.h```: Apply weighting(Trigger, L1Prefiring, Pileup, top pT, ...)
- ```EventSelections.h```: Define signal or control regions
- ```include/ScaleFactors.h```: Calculate scale factors for each event and apply them to each region
- ```include/PlottingBase.h```: Define and fill histograms
