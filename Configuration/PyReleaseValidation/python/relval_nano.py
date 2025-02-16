from  Configuration.PyReleaseValidation.relval_steps import *
import math

workflows = Matrix()

class WFN:
    def __init__(self):
        self.index=0
        self.subindex=0
    def __call__(self):
        offset=2500
        if self.subindex==100:
            print("this is not going to work nicely")
            self.subindex=0/0
        r=float(f'{offset}.{self.index}{self.subindex:02d}')
        self.subindex+=1
        return r
    def next(self):
        self.index+=1
        self.subindex=0
    def subnext(self):
        # go to the next tenth for the subindex 10 because of 02d formating
        self.subindex = math.ceil(self.subindex/10.)*10

_runOnly20events={'-n':'20'}
_run10kevents={'-n':'10000'}

_NANO_data = merge([{'-s':'NANO,DQM:@nanoAODDQM',
                     '--process':'NANO',
                     '--data':'',
                     '--eventcontent':'NANOAOD,DQM',
                     '--datatier':'NANOAOD,DQMIO',
                     '-n':'10000',
                     '--customise':'"Configuration/DataProcessing/Utils.addMonitoring"'
                 }])
_HARVEST_nano = merge([{'-s':'HARVESTING:@nanoAODDQM',
                        '--filetype':'DQM',
                        '--filein':'file:step2_inDQM.root',
                        '--conditions':'auto:run2_data'## this is fake for harvesting
                    }])
_HARVEST_data = merge([_HARVEST_nano, {'--data':''}])


run2_lumis={ 277168: [[1, 1708]],
             277194: [[913, 913], [916, 916], [919, 919], [932, 932], [939, 939]],
             283877: [[1, 1496]],
             299649: [[155, 332]],
             303885: [[60, 2052]],
             305040: [[200, 700]],
             305050: [[200, 700]],
             305234: [[1, 200]],
             305377: [[1, 500]],
             315489: [[1, 100]],
             320822: [[1, 200]],
         }
run3_lumis={}

_NANO_mc = merge([{'-s':'NANO,DQM:@nanoAODDQM',
                   '--process':'NANO',
                   '--mc':'',
                   '--eventcontent':'NANOAODSIM,DQM',
                   '--datatier':'NANOAODSIM,DQMIO',
                   '-n':'10000',
                   '--customise':'"Configuration/DataProcessing/Utils.addMonitoring"'
               }])
_HARVEST_mc = merge([_HARVEST_nano, {'--mc':''}])
steps['HRV_NANO_mc'] = _HARVEST_mc
steps['HRV_NANO_data'] = _HARVEST_data

##8.0 INPUT and workflows
steps['TTbarMINIAOD8.0'] = {'INPUT':InputInfo(location='STD',
                                              dataSet='/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')}
steps['NANO_mc8.0']=merge([{'--era':'Run2_2016,run2_miniAOD_80XLegacy',
                            '--conditions':'auto:run2_mc'},
                           _NANO_mc])

steps['MuonEG2016HMINIAOD8.0'] = {'INPUT':InputInfo(location='STD',ls=run2_lumis,
                                                   dataSet='/MuonEG/Run2016H-03Feb2017_ver2-v1/MINIAOD ')}
steps['NANO_data8.0']=merge([{'--era':'Run2_2016,run2_miniAOD_80XLegacy',
                              '--conditions':'auto:run2_data'},
                             _NANO_data])


##9.4 INPUT and workflows
## de facto, 9.4v1 mini is deprecated
steps['TTbar2016MINIAOD9.4v2'] = {'INPUT':InputInfo(location='STD',
                                                    dataSet='/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM')}
steps['NANO_2016mc9.4v2']=merge([{'--era':'Run2_2016,run2_nanoAOD_94X2016',
                                  '--conditions':'auto:run2_mc'},
                                 _NANO_mc])
steps['TTbarMINIAOD9.4v2'] = {'INPUT':InputInfo(location='STD',
                                                dataSet='/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM')}
steps['NANO_mc9.4v2']=merge([{'--era':'Run2_2017,run2_nanoAOD_94XMiniAODv2',
                              '--conditions':'auto:phase1_2017_realistic'},
                             _NANO_mc])

##17Jul2018 campaign is CMSSW_9_4_9
steps['MuonEG2016MINIAOD9.4v2'] = {'INPUT':InputInfo(location='STD',ls=run2_lumis,
                                                    dataSet='/MuonEG/Run2016H-17Jul2018-v1/MINIAOD')}
steps['NANO_2016data9.4v2']=merge([{'--era':'Run2_2016,run2_nanoAOD_94X2016',
                                    '--conditions':'auto:run2_data'},
                                   _NANO_data])

##17Nov2017 campaign is CMSSW_9_4_0
## de facto, 9.4v1 mini is deprecated

##31Mar2018 campaign is CMSSW_9_4_5_cand1 
steps['MuonEG2017MINIAOD9.4v2'] = {'INPUT':InputInfo(location='STD',ls=run2_lumis,
                                                     dataSet='/MuonEG/Run2017E-31Mar2018-v1/MINIAOD')}
steps['NANO_data9.4v2']=merge([{'--era':'Run2_2017,run2_nanoAOD_94XMiniAODv2',
                                '--conditions':'auto:run2_data'},
                               _NANO_data])

##10.2 INPUT and worflows
steps['TTbarMINIAOD10.2'] = {'INPUT':InputInfo(location='STD',
                                               dataSet='/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM')}
steps['NANO_mc10.2']=merge([{'--era':'Run2_2018,run2_nanoAOD_102Xv1',
                             '--conditions':'auto:phase1_2018_realistic'},
                            _NANO_mc])

##17Sep2018 campaign is CMSSW_10_2_4_patch1 
steps['MuonEG2018MINIAOD10.2'] = {'INPUT':InputInfo(location='STD',ls=run2_lumis,
                                                    dataSet='/MuonEG/Run2018A-17Sep2018-v1/MINIAOD')}
steps['NANO_data10.2']=merge([{'--era':'Run2_2018,run2_nanoAOD_102Xv1',
                               '--conditions':'auto:run2_data'},
                              _NANO_data])

##10.6 INPUT and worflows
steps['TTbarMINIAOD10.6_UL16v1'] = {'INPUT':InputInfo(location='STD',
                                                      dataSet='/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1/MINIAODSIM')}
steps['NANO_mc10.6ul16v1']=merge([{'--era':'Run2_2016,run2_nanoAOD_106Xv1',
                                   '--conditions':'auto:run2_mc'},
                                  _NANO_mc])
steps['TTbarMINIAOD10.6_UL16v2'] = {'INPUT':InputInfo(location='STD',
                                                      dataSet='/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM')}
steps['NANO_mc10.6ul16v2']=merge([{'--era':'Run2_2016,run2_nanoAOD_106Xv2',
                                   '--conditions':'auto:run2_mc'},
                                  _NANO_mc])
##2017 looking Monte-Carlo: two versions in 10.6
steps['TTbarMINIAOD10.6_UL17v1'] = {'INPUT':InputInfo(location='STD',
                                                      dataSet='/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v2/MINIAODSIM')}
steps['NANO_mc10.6ul17v1']=merge([{'--era':'Run2_2017,run2_nanoAOD_106Xv1',
                                   '--conditions':'auto:phase1_2017_realistic'},
                                  _NANO_mc])
steps['TTbarMINIAOD10.6_UL17v2'] = {'INPUT':InputInfo(location='STD',
                                                      dataSet='/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM')}
steps['NANO_mc10.6ul17v2']=merge([{'--era':'Run2_2017,run2_nanoAOD_106Xv2',
                                   '--conditions':'auto:phase1_2017_realistic'},
                                  _NANO_mc])

steps['TTbarMINIAOD10.6_UL18v1'] = {'INPUT':InputInfo(location='STD',
                                                      dataSet='/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1/MINIAODSIM')}
steps['NANO_mc10.6ul18v1']=merge([{'--era':'Run2_2018,run2_nanoAOD_106Xv1',
                                   '--conditions':'auto:phase1_2018_realistic'},
                                  _NANO_mc])
steps['TTbarMINIAOD10.6_UL18v2'] = {'INPUT':InputInfo(location='STD',
                                                      dataSet='/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM')}
steps['NANO_mc10.6ul18v2']=merge([{'--era':'Run2_2018,run2_nanoAOD_106Xv2',
                                   '--conditions':'auto:phase1_2018_realistic'},
                                  _NANO_mc])

##21Feb2020 campaign is CMSSW_10_6_8_patch1 
steps['MuonEG2016MINIAOD10.6v1'] = {'INPUT':InputInfo(location='STD',ls=run2_lumis,
                                                      dataSet='/MuonEG/Run2016E-21Feb2020_UL2016_HIPM-v1/MINIAOD')}
steps['NANO_data10.6ul16v1']=merge([{'--era':'Run2_2016,run2_nanoAOD_106Xv1,tracker_apv_vfp30_2016',
                                     '--conditions':'auto:run2_data'},
                                    _NANO_data])
##09Aug2019 campaign is CMSSW_10_6_2 
steps['MuonEG2017MINIAOD10.6v1'] = {'INPUT':InputInfo(location='STD',ls=run2_lumis,
                                                      dataSet='/MuonEG/Run2017F-09Aug2019_UL2017-v1/MINIAOD')}
steps['NANO_data10.6ul17v1']=merge([{'--era':'Run2_2017,run2_nanoAOD_106Xv1',
                                     '--conditions':'auto:run2_data'},
                                    _NANO_data])

##12Nov2019 campaign is CMSSW_10_6_4_patch1 
steps['MuonEG2018MINIAOD10.6v1'] = {'INPUT':InputInfo(location='STD',ls=run2_lumis,
                                                      dataSet='/MuonEG/Run2018A-12Nov2019_UL2018_rsb-v1/MINIAOD')}
steps['NANO_data10.6ul18v1']=merge([{'--era':'Run2_2018,run2_nanoAOD_106Xv1',
                                     '--conditions':'auto:run2_data'},
                                    _NANO_data])

##HIPM_UL2016_MiniAODv2 campaign is CMSSW_10_6_25
steps['MuonEG2016MINIAOD10.6v2'] = {'INPUT':InputInfo(location='STD',ls=run2_lumis,
                                                      dataSet='/MuonEG/Run2016E-HIPM_UL2016_MiniAODv2-v2/MINIAOD')}
steps['NANO_data10.6ul16v2']=merge([{'--era':'Run2_2016,run2_nanoAOD_106Xv2,tracker_apv_vfp30_2016',
                                     '--conditions':'auto:run2_data'},
                                    _NANO_data])
##UL2017_MiniAODv2 campaign is CMSSW_10_6_20 
steps['MuonEG2017MINIAOD10.6v2'] = {'INPUT':InputInfo(location='STD',ls=run2_lumis,
                                                      dataSet='/MuonEG/Run2017F-UL2017_MiniAODv2-v1/MINIAOD')}
steps['NANO_data10.6ul17v2']=merge([{'--era':'Run2_2017,run2_nanoAOD_106Xv2',
                                     '--conditions':'auto:run2_data'},
                                    _NANO_data])

##UL2018_MiniAODv2 campaign is CMSSW_10_6_20 
steps['MuonEG2018MINIAOD10.6v2'] = {'INPUT':InputInfo(location='STD',ls=run2_lumis,
                                                      dataSet='/MuonEG/Run2018D-UL2018_MiniAODv2-v1/MINIAOD')}
steps['NANO_data10.6ul18v2']=merge([{'--era':'Run2_2018,run2_nanoAOD_106Xv2',
                                     '--conditions':'auto:run2_data'},
                                    _NANO_data])
##12.2 INPUT (mc only)
steps['TTbarMINIAOD12.2'] = {'INPUT':InputInfo(location='STD',
                                               dataSet='/TTToSemiLeptonic_TuneCP5_13p6TeV-powheg-pythia8/Run3Winter22MiniAOD-FlatPU0to70_122X_mcRun3_2021_realistic_v9-v2/MINIAODSIM')}
steps['NANO_mc12.2_v10']=merge([{'--era':'Run3,run3_nanoAOD_122',
                                 '--conditions':'auto:phase1_2022_realistic',
                                 '-s':'NANO:PhysicsTools/NanoAOD/V10/nano_cff,DQM:@nanoAODDQM'},
                                _NANO_mc])
steps['NANO_mc12.2']=merge([{'--era':'Run3,run3_nanoAOD_122',
                             '--conditions':'auto:phase1_2022_realistic'},
                            _NANO_mc])

##12.4 INPUT 
steps['TTbarMINIAOD12.4'] = {'INPUT':InputInfo(location='STD',
                                               ## to be updated as soon as some TTbar appears in a 12.4 campaign
                                               dataSet='/MinBias_TuneCP5_14TeV-pythia8/Run3Summer22MiniAODv3-NoPU_Pilot_124X_mcRun3_2022_realistic_v11-v2/MINIAODSIM ')}
steps['NANO_mc12.4_v10']=merge([{'--era':'Run3,run3_nanoAOD_124',
                                 '--conditions':'auto:phase1_2022_realistic',
                                 '-s':'NANO:PhysicsTools/NanoAOD/V10/nano_cff,DQM:@nanoAODDQM'},
                                _NANO_mc])
steps['NANO_mc12.4']=merge([{'--era':'Run3,run3_nanoAOD_124',
                             '--conditions':'auto:phase1_2022_realistic'},
                            _NANO_mc])

steps['MuonEG2022MINIAOD12.4'] = {'INPUT':InputInfo(location='STD',ls=run3_lumis,
                                                    dataSet='/MuonEG/Run2022D-PromptReco-v2/MINIAOD')}
steps['NANO_data12.4_v10']=merge([{'--era':'Run3,run3_nanoAOD_124',
                                   '--conditions':'auto:run3_data',
                                   '-s':'NANO:PhysicsTools/NanoAOD/V10/nano_cff,DQM:@nanoAODDQM'},
                                  _NANO_data])
steps['NANO_data12.4']=merge([{'--era':'Run3,run3_nanoAOD_124',
                               '--conditions':'auto:run3_data'},
                              _NANO_data])

##12.6 workflows ("from scratch")
steps['TTBarMINIAOD12.6'] = {'INPUT':InputInfo(location='STD',ls=run3_lumis,
                                               ## this is a dataset from the last pre-release: to be updated much too often IMO
                                               dataSet='/RelValTTbar_14TeV/CMSSW_12_6_0_pre2-PU_125X_mcRun3_2022_realistic_v3-v1/MINIAODSIM')}
steps['NANO_mc12.6_v10']=merge([{'--era':'Run3',
                                 '--conditions':'auto:phase1_2022_realistic',
                                 '-s':'NANO:PhysicsTools/NanoAOD/V10/nano_cff,DQM:@nanoAODDQM'},
                                _NANO_data])
steps['NANO_mc12.6']=merge([{'--era':'Run3',
                             '--conditions':'auto:phase1_2022_realistic'},
                            _NANO_data])


_wfn=WFN()
################
#8.X input
workflows[_wfn()] = ['mc80X', ['TTbarMINIAOD8.0','NANO_mc8.0', 'HRV_NANO_mc']]
workflows[_wfn()] = ['data80X', ['MuonEG2016HMINIAOD8.0','NANO_data8.0', 'HRV_NANO_data']]

_wfn.next()
################
#9.4 input
#workflows[_wfn()] = ['mc94X', ['TTbarMINIAOD9.4v1','NANO_mc9.4v1']] ##deprecated input
workflows[_wfn()] = ['mc94X2016', ['TTbar2016MINIAOD9.4v2','NANO_2016mc9.4v2', 'HRV_NANO_mc']]
workflows[_wfn()] = ['mc94Xv2', ['TTbarMINIAOD9.4v2','NANO_mc9.4v2', 'HRV_NANO_mc']]
_wfn.subnext()
#workflows[_wfn()] = ['data94X', ['MuonEG2017MINIAOD9.4v1','NANO_data9.4v1']] ## deprecated
_wfn.subnext()
workflows[_wfn()] = ['data94X2016', ['MuonEG2016MINIAOD9.4v2','NANO_2016data9.4v2', 'HRV_NANO_data']]
_wfn.subnext()
workflows[_wfn()] = ['data94Xv2', ['MuonEG2017MINIAOD9.4v2','NANO_data9.4v2', 'HRV_NANO_data']]

_wfn.next()
################
#10.2 input
workflows[_wfn()] = ['mc102X', ['TTbarMINIAOD10.2','NANO_mc10.2', 'HRV_NANO_mc']]
_wfn.subnext()
workflows[_wfn()] = ['data102X', ['MuonEG2018MINIAOD10.2','NANO_data10.2', 'HRV_NANO_data']]

_wfn.next()
################
#10.6 input
workflows[_wfn()] = ['mc106Xul16', ['TTbarMINIAOD10.6_UL16v1','NANO_mc10.6ul16v1', 'HRV_NANO_mc']]
workflows[_wfn()] = ['mc106Xul17', ['TTbarMINIAOD10.6_UL17v1','NANO_mc10.6ul17v1', 'HRV_NANO_mc']]
workflows[_wfn()] = ['mc106Xul18', ['TTbarMINIAOD10.6_UL18v1','NANO_mc10.6ul18v1', 'HRV_NANO_mc']]
_wfn.subnext()
workflows[_wfn()] = ['mc106Xul17v2', ['TTbarMINIAOD10.6_UL16v2','NANO_mc10.6ul16v2', 'HRV_NANO_mc']]
workflows[_wfn()] = ['mc106Xul17v2', ['TTbarMINIAOD10.6_UL17v2','NANO_mc10.6ul17v2', 'HRV_NANO_mc']]
workflows[_wfn()] = ['mc106Xul17v2', ['TTbarMINIAOD10.6_UL18v2','NANO_mc10.6ul18v2', 'HRV_NANO_mc']]
_wfn.subnext()
workflows[_wfn()] = ['data106Xul16', ['MuonEG2016MINIAOD10.6v1', 'NANO_data10.6ul16v1', 'HRV_NANO_data']]
workflows[_wfn()] = ['data106Xul17', ['MuonEG2017MINIAOD10.6v1', 'NANO_data10.6ul17v1', 'HRV_NANO_data']]
workflows[_wfn()] = ['data106Xul18', ['MuonEG2018MINIAOD10.6v1', 'NANO_data10.6ul18v1', 'HRV_NANO_data']]
_wfn.subnext()
workflows[_wfn()] = ['data106Xul16v2', ['MuonEG2016MINIAOD10.6v2', 'NANO_data10.6ul16v2', 'HRV_NANO_data']]
workflows[_wfn()] = ['data106Xul17v2', ['MuonEG2017MINIAOD10.6v2', 'NANO_data10.6ul17v2', 'HRV_NANO_data']]
workflows[_wfn()] = ['data106Xul18v2', ['MuonEG2018MINIAOD10.6v2', 'NANO_data10.6ul18v2', 'HRV_NANO_data']]

_wfn.next()
################
#12.2 input
workflows[_wfn()] = ['mc122Xrun3_v10', ['TTbarMINIAOD12.2','NANO_mc12.2_v10', 'HRV_NANO_mc']]
workflows[_wfn()] = ['mc122Xrun3', ['TTbarMINIAOD12.2','NANO_mc12.2', 'HRV_NANO_mc']]

_wfn.next()
################
#12.4 input
## these are borken because of tau configuration in NANO ATM: they should be re-enabled when a fix gets in
workflows[_wfn()] = ['mc124Xrun3_v10', ['TTbarMINIAOD12.4','NANO_mc12.4_v10', 'HRV_NANO_mc']]
workflows[_wfn()] = ['mc124Xrun3', ['TTbarMINIAOD12.4','NANO_mc12.4', 'HRV_NANO_mc']]
_wfn.subnext()
workflows[_wfn()] = ['data124Xrun3_v10', ['MuonEG2022MINIAOD12.4','NANO_data12.4_v10', 'HRV_NANO_data']]
workflows[_wfn()] = ['data124Xrun3', ['MuonEG2022MINIAOD12.4','NANO_data12.4', 'HRV_NANO_data']]

_wfn.next()
################
#12.6 workflows
## these two workflows should be creating a sample "from scratch" instead of using a pre-release sample as input
workflows[_wfn()] = ['mc126X_v10', ['TTBarMINIAOD12.6','NANO_mc12.6_v10', 'HRV_NANO_mc']]
workflows[_wfn()] = ['mc126X', ['TTBarMINIAOD12.6','NANO_mc12.6', 'HRV_NANO_mc']]
_wfn.subnext()

_wfn.next()
################
