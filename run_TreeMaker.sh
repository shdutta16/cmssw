#!/bin/bash


cmsRun EDAnalyzers/TreeMaker/python/ConfFile_cfg.py \
    sourceFile=\
sourceFiles/TTToSemiLepton_TuneCP5_14TeV-powheg-pythia8_Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2_GEN-SIM-DIGI-RAW/TTToSemiLepton_TuneCP5_14TeV-powheg-pythia8_Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2_GEN-SIM-DIGI-RAW.txt \
    genEleFilter=0 \
    genPartonFilter=1 \
    isGunSample=0 \
    modTICLele=1 \
    modTICLeleWithRerunTICL=1 \
    rerunTICL=1 \
    onRaw=1 \
    maxEvents=10 \
