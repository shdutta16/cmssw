#!/bin/bash


cmsRun EDAnalyzers/TreeMaker/python/ConfFile_cfg.py \
    sourceFile=\
sourceFiles/QCD_Pt-15to3000_TuneCP5_Flat_14TeV-pythia8_Phase2HLTTDRWinter20DIGI-PU200_castor_110X_mcRun4_realistic_v3-v2_GEN-SIM-DIGI-RAW/QCD_Pt-15to3000_TuneCP5_Flat_14TeV-pythia8_Phase2HLTTDRWinter20DIGI-PU200_castor_110X_mcRun4_realistic_v3-v2_GEN-SIM-DIGI-RAW.txt \
    genEleFilter=0 \
    genPartonFilter=1 \
    isGunSample=0 \
    modTICLele=1 \
    modTICLeleWithRerunTICL=1 \
    rerunTICL=1 \
    onRaw=1 \
    printTime=1 \
    eventRange="1:4538-1:4538" \
    maxEvents=1 \
