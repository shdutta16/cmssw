#!/bin/bash

cmsRun EDAnalyzers/TreeMaker/python/ConfFile_cfg.py \
    sourceFile=\
sourceFiles/QCD_Pt-50to80_EMEnriched_TuneCP5_14TeV_pythia8_Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2_GEN-SIM-DIGI-RAW/QCD_Pt-50to80_EMEnriched_TuneCP5_14TeV_pythia8_Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2_GEN-SIM-DIGI-RAW.txt \
    genEleFilter=0 \
    genPartonFilter=0 \
    isGunSample=0 \
    modTICLele=0 \
    modTICLeleWithRerunTICL=0 \
    rerunTICL=0 \
    onRaw=1 \
    debugFile=0 \
    maxEvents=50 \
