<pre>


1. Running the tree maker on RAW:
cmsRun EDAnalyzers/TreeMaker/python/ConfFile_cfg.py \
    onRaw=1 \
    inputFiles=[INPUTFILES] \
    sourceFile=[file containing list of input files. Will use only if <inputFiles> is not provided] \
    modTICLele=1 \
    rerunTICL=[1 if you want to rerun the TICL sequence] \
    modTICLeleWithRerunTICL=[1 if you want TICL-electrons with the rerun TICL objects]




##############################
Git stuff
##############################

1. Push changes:
    a) git add .
    b) git commit
    c) git push -u my-cmssw SohamBhattacharya/HGCal_ele-reco:HGCal_ele-reco 

2. Delete from remote branch ONLY:
    git rm --cached <file>
    git rm -r --cached <directory>
    <1b>
    <1c>



</pre>
