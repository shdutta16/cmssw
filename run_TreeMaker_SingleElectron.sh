#!/bin/bash


cmsRun EDAnalyzers/TreeMaker/python/ConfFile_cfg.py \
    sourceFile=\
sourceFiles/SingleElectron_PT2to200_Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3_ext2-v2_GEN-SIM-modRECO/SingleElectron_PT2to200_Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3_ext2-v2_GEN-SIM-modRECO.txt \
    genEleFilter=1 \
    genPartonFilter=0 \
    isGunSample=1 \
    modTICLele=1 \
    modTICLeleWithRerunTICL=1 \
    rerunTICL=1 \
    onRaw=0 \
    maxEvents=10 \
