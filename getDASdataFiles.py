import os


outDir = "sourceFiles"

prefix = "root://cms-xrd-global.cern.ch//store"
toReplace = "/store"


l_sampleName = [
    
    #"/QCD_Pt-15to3000_TuneCP5_Flat_14TeV-pythia8/Phase2HLTTDRWinter20DIGI-PU200_castor_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW",
    
    #"/SingleElectron_PT2to200/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext2-v1/FEVT",
    #"/SingleElectron_PT200to500/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext2-v1/FEVT",
    "/DoubleElectron_FlatPt-1To100/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v2/GEN-SIM-DIGI-RAW-MINIAOD",
    
    #"/SinglePhoton_PT2to200/Phase2HLTTDRWinter20DIGI-NoPU_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW",
    #"/SinglePhoton_PT2to200/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext2-v1/FEVT",
    "/SinglePhoton_PT200to500/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext2-v1/FEVT",
    "/DoublePhoton_FlatPt-1To100/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext1-v2/FEVT",
    
    #"/TTToSemiLepton_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/FEVT",
]


for iSample in range(0, len(l_sampleName)) :
    
    sampleName = l_sampleName[iSample]
    
    sampleName_mod = sampleName[1:].replace("/", "_")
    
    outDir_mod = "%s/%s" %(outDir, sampleName_mod)
    
    command = "mkdir -p %s" %(outDir_mod)
    print "Command:", command
    print ""
    os.system(command)
    
    outFile = "%s/%s.txt" %(outDir_mod, sampleName_mod)
    
    command = "dasgoclient -query=\"file dataset=%s\" > %s" %(sampleName, outFile)
    print "Command:", command
    print ""
    os.system(command)
    
    fileContent = ""
    
    print "Replacing \"%s\" with \"%s\" in file." %(toReplace, prefix)
    print ""
    
    with open(outFile, "r") as f :
        
        fileContent = f.read()
    
    fileContent = fileContent.replace(toReplace, prefix)
    
    with open(outFile, "w") as f :
        
        f.write(fileContent)
