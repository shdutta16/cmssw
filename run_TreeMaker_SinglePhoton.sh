#!/bin/bash


cmsRun EDAnalyzers/TreeMaker/python/ConfFile_cfg.py \
    sourceFile=\
sourceFiles/SinglePhoton_PT2to200_Phase2HLTTDRWinter20DIGI-NoPU_110X_mcRun4_realistic_v3-v2_GEN-SIM-DIGI-RAW/SinglePhoton_PT2to200_Phase2HLTTDRWinter20DIGI-NoPU_110X_mcRun4_realistic_v3-v2_GEN-SIM-DIGI-RAW.txt \
    genEleFilter=0 \
    genPhoFilter=1 \
    genPartonFilter=0 \
    isGunSample=1 \
    modTICLele=0 \
    modTICLeleWithRerunTICL=0 \
    rerunTICL=0 \
    onRaw=1 \
    debugFile=0 \
    maxEvents=2000 \
    outputDir="output/SinglePhoton_PT2to200_Phase2HLTTDRWinter20DIGI-NoPU_110X_mcRun4_realistic_v3-v2_GEN-SIM-DIGI-RAW" \


#cmsRun EDAnalyzers/TreeMaker/python/ConfFile_cfg.py \
#    sourceFile=\
#sourceFiles/SinglePhoton_PT2to200_Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext2-v1_FEVT/SinglePhoton_PT2to200_Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext2-v1_FEVT.txt \
#    genEleFilter=0 \
#    genPhoFilter=1 \
#    genPartonFilter=0 \
#    isGunSample=1 \
#    modTICLele=0 \
#    modTICLeleWithRerunTICL=0 \
#    rerunTICL=0 \
#    onRaw=1 \
#    debugFile=0 \
#    maxEvents=10 \
