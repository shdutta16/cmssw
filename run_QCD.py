from __future__ import print_function

import argparse
import numpy
import os


cwd = os.getcwd()


# Argument parser
parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter)

# TreeMaker_SingleElectron_PT2to100_PhaseIITDRSpring19DR-PU200_106X_upgrade2023_realistic_v3-v1_GEN-SIM-DIGI-RAW
parser.add_argument(
    "--samples",
    help = "Names of the samples",
    type = str,
    nargs = "*",
    required = True,
)

parser.add_argument(
    "--submit",
    help = "Submit jobs for the samples",
    action = "store_true",
    default = False,
)


# Parse arguments
args = parser.parse_args()


for iSamp, sampleName in enumerate(args.samples) :
    
    print("*"*100)
    print("Sample %d/%d: %s" %(iSamp+1, len(args.samples), sampleName))
    print("*"*100)
    
    processName = "TreeMaker_%s" %(sampleName)
    sourceFileName = "sourceFiles/%s/%s.txt" %(sampleName, sampleName)
    
    cmd_str = (
        "python run_condor.py "
        "--processName {processName} "
        "--inputFileList {sourceFileName} "
        "--cmsRunFile EDAnalyzers/TreeMaker/python/ConfFile_cfg.py "
        "--outputDir /eos/cms/store/group/phys_egamma/sobhatta/HGCal_TreeMaker "
        "--nUnitPerJob 1 "
        "--cmsRunOptions \"genPartonFilter=1 isGunSample=0 onRaw=1\" "
        "--suffix \"\" "
    ).format(
        processName = processName,
        sourceFileName = sourceFileName,
    )
    
    print("Command:")
    print(cmd_str)
    
    if (args.submit) :
        
        cmd_retVal = 99
        
        print("Submitting...")
        
        while (cmd_retVal) :
            
            cmd_retVal = os.system(cmd_str)
        
    print("\n")
