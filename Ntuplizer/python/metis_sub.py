import time

from metis.Sample import DBSSample,DirectorySample
from metis.CMSSWTask import CMSSWTask
from metis.CondorTask import CondorTask

from metis.StatsParser import StatsParser

if __name__ == "__main__":

    dataset_names = [
            "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM",
            # "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM",

            # "/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
            # "/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v3/MINIAODSIM",
            # "/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
            # "/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
            # "/WWTo2L2Nu_DoubleScattering_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
            # "/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
            # "/ZZTo4L_13TeV_powheg_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",

            # "/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
            # "/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",

            # "/DoubleEG/Run2016H-03Feb2017_ver2-v1/MINIAOD",
            # "/DoubleMuon/Run2016H-03Feb2017_ver2-v1/MINIAOD",
            # "/MuonEG/Run2016H-03Feb2017_ver2-v1/MINIAOD",

            ]

    cmssw_tasks = []
    for dsname in dataset_names:
        cmsswver = "CMSSW_8_0_28"
        pset = "run_cfg.py"
        # tarfile = "tarfile.tar.gz"
        # tarfile = "tarfile_15x29.tar.gz"
        # tarfile = "tarfile_15x29_tracksv2.tar.gz"
        # tarfile = "tarfile_15x29_timev1.tar.gz"
        tarfile = "tarfile_15x29_sconlyv1.tar.gz"
        # tag = "v12"
        # tag = "v13" # 15x29
        # tag = "v14" # 29x29
        # tag = "v15" # back to 15x29 ttbar, dy w/ track
        # tag = "v16" # add R to track variables
        # tag = "v17" # add rhs_time, rhs_chi2 branches, DY only
        tag = "v18" # DY only, SC only instead of all towers
        pset_args = "print"
        if "Run2016" in dsname:
            pset_args = "data=True"
        cmssw_task = CMSSWTask(
                sample = DBSSample( dataset=dsname ),
                executable = "condor_executable.sh",
                special_dir = "elemva/ProjectMetis",
                open_dataset = False,
                events_per_output = 700000,
                output_name = "output.root",
                output_is_tree = False,
                check_expectedevents = True,
                tag = tag,
                pset = pset,
                cmssw_version = cmsswver,
                tarfile = tarfile,
                pset_args = pset_args,
                publish_to_dis = False,
                )
        cmssw_tasks.append(cmssw_task)

    for i in range(100):

        total_summary = {}
            
        for cmssw_task in cmssw_tasks:
            cmssw_task.process()

            total_summary[cmssw_task.get_sample().get_datasetname()] = cmssw_task.get_task_summary()

        # parse the total summary and write out the dashboard
        StatsParser(data=total_summary, webdir="~/public_html/dump/metis/").do()

        # 1 hr power nap
        time.sleep(10.*60)

