<pre>


1. Running the tree maker on RAW:
cmsRun EDAnalyzers/TreeMaker/python/ConfFile_cfg.py \
    onRaw=1 \
    inputFiles=[INPUTFILES] \
    sourceFile=[file containing list of input files. Will use only if <inputFiles> is not provided] \
    modTICLele=1 \
    rerunTICL=[1 if you want to rerun the TICL sequence] \
    modTICLeleWithRerunTICL=[1 if you want TICL-electrons with the rerun TICL objects] \
    TICLeleGenMatchDR=[TICLELEGENMATCHDR (by default all TICL-electrons will be stored as the default value is very large)]




##############################
Git stuff
##############################

1. Push changes:
    a) git add .
    b) git commit
    c) git push -u my-cmssw from-CMSSW_11_1_0_pre3:HGCal_ele-reco_analysis_11_1_0_pre3 

2. Delete from remote branch ONLY:
    git rm --cached <file>
    git rm -r --cached <directory>
    <1b>
    <1c>

3. If files become untracked for some weird reason, do:
    for f in $(find | grep -v .git | grep -v .pyc | grep -v __init__); do echo $f; git update-index --no-skip-worktree $f; done

</pre>
