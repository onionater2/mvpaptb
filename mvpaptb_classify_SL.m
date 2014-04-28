function [subj, printregressor]=mvpaptb_classify_SL(study, task, runsincluded, mvparootdir, saveimgsdir, subjectID, subj, binarized, convolved, condnames, mask_description, xvalselector, averagingselector, averagingselectorcomb, other_args, class_args)
% created by AES 4/18/13
% called by run_classification_EIB
% takes as input a single binary regressor and convolved regressor (need to
% run separate times for separate classifications)
%


%% other fixed parameters
rootdir=['/mindhive/saxelab2/'];
studydir=[rootdir study '/'];
% Check to make sure the Neuralnetwork toolbox is in the path or this
% won't work.
%TO DO: update to check for toolboxes for other analyses
if ~exist('newff') %#ok<EXIST>
    error('This tutorial requires the neural networking toolbox, if it is unavailable this will not execute');
end


subjdir=[studydir subjectID '/'];
cd(subjdir)

load([mvparootdir 'subjinfo_' task '_' runsincluded]); %% contains variables maskname, mvpa_mask_description, boldnames, condregressors, condnames, runselectors
load([mvparootdir 'discriminations_' task '_' runsincluded]); %contains discriminations

condregressors_binarized=binarized;
condregressors_convolved=convolved;
blocklabels=eval(averagingselector);
blocklabelscomb=eval(averagingselectorcomb);


% initialize and set the selectors objectt
%% run based selectors (don't need shifting)
subj = init_object(subj,'selector','runs');
subj = set_mat(subj,'selector','runs',runselectors);

subj = init_object(subj,'selector','evenodd');
subj = set_mat(subj,'selector','evenodd',evenoddselectors);


% initialize and set contents of the regressors object in the subj structure,fit and add a
% cell array of condnames to the object for future reference
% these regressors are based on onsets/durations from SPM.Sess.U

regressorlabel='conds_binarized';
subj = init_object(subj,'regressors',regressorlabel);
subj = set_mat(subj,'regressors',regressorlabel,eval(condregressors_binarized));
subj = set_objfield(subj,'regressors',regressorlabel,'condnames',eval(condnames));


% shift regressor by hemodynamic delay (in TR units) 
shift=other_args.hemodynamic_delay/TR;
subj = shift_regressors(subj,regressorlabel,'runs', shift);

% add convolved regressor, which was calculated by convolving HRF to the above regressors (of course, don't need to shift this since convolution does
% that). currently just using binary for everything though.
subj = init_object(subj,'regressors','conds_convolved');
subj = set_mat(subj,'regressors','conds_convolved',eval(condregressors_convolved));
subj = set_objfield(subj,'regressors','conds_convolved','condnames',eval(condnames));

% initialize and assign additional selectors
% shift each of these by hemodynamic_delay/TR so that it aligns with the shifted binarized
% regressor

artselectors=[zeros(shift,1); artselectors(1:end-shift)];
subj = init_object(subj,'selector','artifacts');
subj = set_mat(subj,'selector','artifacts',artselectors);

restselectors=[zeros(shift,1); restselectors(1:end-shift)];
subj = init_object(subj,'selector','rests');
subj = set_mat(subj,'selector','rests',restselectors);

combselectors=[zeros(shift,1); combselectors(1:end-shift)];
subj = init_object(subj,'selector','combined_art_and_rest');
subj = set_mat(subj,'selector','combined_art_and_rest',combselectors);

%if you are using an xval selector other than runs or evenodd (e.g. one
%based on the events within a run), shift that too
if ~strcmp(xvalselector,'runs') && ~strcmp(xvalselector,'evenodd')
    xvalsel=eval(xvalselector);
    xvalsel=[zeros(shift,1); xvalsel(1:end-shift)];
    subj = init_object(subj,'selector',xvalselector);
    subj = set_mat(subj,'selector', xvalselector,xvalsel); 
end

%% these blocklabels selectors assign different value to each trial. used to determine averaging to create examples
blocklabels=[zeros(shift,1); blocklabels(1:end-shift)];
subj = init_object(subj,'selector','blocklabels');
subj = set_mat(subj,'selector','blocklabels',blocklabels);

%% this selector drops any trials for which a single TR is excluded in art
blocklabelscomb=[zeros(shift,1); blocklabelscomb(1:end-shift)];
subj = init_object(subj,'selector','blocklabelscomb_art_and_rest');
subj = set_mat(subj,'selector','blocklabelscomb_art_and_rest',blocklabelscomb);

timecourseselector=runselectors*0+1; %make vector of ones of length of timecourse
subj = init_object(subj,'selector','wholetimecourse');
subj = set_mat(subj,'selector','wholetimecourse',timecourseselector);

%option to detrend across whole experiment
if other_args.detrend
subj = detrend_pattern(subj, other_args.bolds, 'wholetimecourse');
other_args.bolds=[other_args.bolds '_dt'];
end
%option to highpass filter runs (could implement these in run_classification?)
if other_args.hpfilter
subj = hpfilter_runs(subj,other_args.bolds,'runs', 128,2);
other_args.bolds=[other_args.bolds '_hp'];
end


%% some silliness to set up some conventionalized naming for classification analyses

if other_args.binary
   shiftedregressor=[regressorlabel '_sh' num2str(shift)];
   printregressor=binarized;
else
   shiftedregressor='conds_convolved'; 
   printregressor=convolved;
end

if other_args.averaged
    if other_args.zscore
    patterns=[other_args.bolds '_z_avg'];
    else if ~other_args.zscore
    patterns=[other_args.bolds '_avg'];
    end
    end
    finalregressor=[shiftedregressor '_avg'];
   if other_args.wart
       averagingselector='blocklabelscomb_art_and_rest';
       useselector=[xvalselector '_avg'];
   else
       averagingselector='blocklabels';
       useselector=[xvalselector '_avg'];
   end
else
    if other_args.zscore
    patterns=[other_args.bolds '_z'];
    else if ~other_args.zscore
    patterns=other_args.bolds;
    end
    end
    finalregressor=shiftedregressor;
    if other_args.wart
       additionalselector='combined_art_and_rest';
       useselector=xvalselector;
    else
       additionalselector='rests';
       useselector=xvalselector;
   end
end


%%subj = load_spm_pattern(subj,bolds,mask_description, boldnames); %% CHANGE: now this step occurs in
%%run_classification_EIB because it is time consuming and we run many
%%classifications per ROI
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRE-PROCESSING - z-scoring in time and no-peeking anova

% we want to z-score the EPI data (called 'epi'),
% individually on each run (using the 'runs' selectors)
if other_args.zscore
subj = zscore_runs(subj,other_args.bolds,'runs'); %note zscoring is occuring over full timecourse, not just examples
    fslabel='_z';
    if other_args.averaged
    disp('creating averaged examples...')
    fslabel='_z_avg'; %so that below you grab the right bolds for feature selection
    subj = average_object(subj,'pattern',[other_args.bolds '_z'],averagingselector);
    subj = average_object(subj,'regressors',shiftedregressor,averagingselector);
    subj = average_object(subj,'selector',xvalselector,averagingselector); %% because art and rest are included in blocklabel, run_avg will now exclude these?
    end
else if ~other_args.zscore
    fslabel='';
    if other_args.averaged
        fslabel='_avg'; %so that below you grab the right bolds for feature selection
    subj = average_object(subj,'pattern',other_args.bolds, averagingselector);
    subj = average_object(subj,'regressors',shiftedregressor,averagingselector);
    subj = average_object(subj,'selector',xvalselector,averagingselector); %% because art and rest are included in blocklabel, run_avg will now exclude these?
    end
    end
end

% now, create selector indices for the n different iterations of
% the nminusone %ignoring jumbled runs to enable cross validation on just
% even/odd, for example
disp('creating xval selector indices...')
if other_args.averaged
    crossvalidation_selector=[useselector '_xval'];
    subj = create_xvalid_indices(subj,useselector, 'ignore_jumbled_runs', 'true');
else
    crossvalidation_selector=[useselector '_' additionalselector '_xval'];
    subj = create_xvalid_indices(subj,useselector, 'actives_selname',additionalselector,'new_selstem', crossvalidation_selector, 'ignore_jumbled_runs', 'true');
end

        % some feature selection functions (e.g. activity) need to know the number of xval folds you'll be using
        % so, find the selectors within the crossvalidiation_selector group
        selnames = find_group(subj,'selector',crossvalidation_selector);
        numfolds = length(selnames);

%get anova or activity maps that you can feed into searchlight analysis
if strcmp(other_args.fsfunc,'anova')
    fspatterns=[other_args.bolds fslabel];
    fsregressor=finalregressor;
    fsselector=crossvalidation_selector;
    if other_args.fstype==1 fixed=1; else fixed=0; end
    subj=feature_select_anova(subj, fspatterns,fsregressor,fsselector, numfolds, fixed, other_args.fixednum, 'thresh',other_args.fsthreshold); %generate discriminability map that you can feed into next searchlight feature selection
    fs_args.initialfs_pat_map= find_group(subj,'pattern',[fspatterns '_anova']);
else if strcmp(other_args.fsfunc,'activity')
        if other_args.zscore
        fspatterns=[other_args.bolds '_z'];
        else
        fspatterns=other_args.bolds;
        end
        fsregressor=finalregressor;
        fsselector='rests';
        if other_args.fstype==1 fixed=1; else fixed=0; end
        subj=feature_select_activity(subj, fspatterns,fsregressor,fsselector, numfolds, fixed, other_args.fixednum, 'thresh',other_args.fsthreshold); %generate activity map
        fs_args.initialfs_pat_map=[fspatterns '_anova'];
    end
end
if other_args.fstype==3 fs_args.initialfs_pat_map='none'; end;
fs_args.fixed=fixed;

% specify and run the searchlight 
class_args.train_funct_name=['train_' other_args.classifier];
class_args.test_funct_name=['test_' other_args.classifier];
class_args.ignore_1ofn = 'false';
if exist('class_args.penalty')
    class_args.penalty=class_args.penalty*floor(other_args.searchlightnumvox);
end
statmap_srch_arg.adj_list = subj.adj_sphere;
statmap_srch_arg.obj_funct = 'statmap_classify_AES';
statmap_srch_arg.scratch.class_args = class_args;
statmap_srch_arg.scratch.perfmet_funct = 'perfmet_maxclass';
statmap_srch_arg.scratch.perfmet_args = struct([]);

fs_args.fstype=other_args.fstype;
fs_args.fixednum=other_args.fixednum;
fs_args.voxelrange=other_args.voxelrange;

[subj, outputpatname] = feature_select_SL(subj,patterns,finalregressor,crossvalidation_selector, numfolds, fs_args, 'statmap_funct','statmap_searchlight_AES', 'statmap_arg', statmap_srch_arg, 'new_map_patname',[patterns '_srch'], 'thresh',other_args.fsthreshold);
[searchlightmaps] = find_group(subj,'pattern',outputpatname);

tempmat=get_mat(subj, 'regressors', finalregressor);
numCondsInDisc=size(tempmat,1); %figure out number of conds based on regressor matrix
chancelevel=1/numCondsInDisc; %compute chance level to subtract voxel-wise

%now lets make a brain image!
patheader=subj.patterns{1}.header; % go back to that first input pattern (assumes the first pat is your raw bold imgs) and grab the header. to do: would be prettier to actually grab based on epi name...
patheader.vol=patheader.vol(1); %that was a 4d pattern but we only need the .vol struct for the first volum
patheader.vol=rmfield(patheader.vol{1}, 'pinfo'); %don't want scaling from the bold
%patheader.vol={patheader.vol};
%label the imgs by the training/testing data
for fold=1:length(searchlightmaps)
   if ~strcmp(xvalselector, 'runs')
       if strcmp(xvalselector,'crossmatchedruns')
      trainset='trainruns';
       else
       for foldsearch=1:length(other_args.foldlabels)
           foldname=other_args.foldlabels(foldsearch).name
           if strcmp(xvalselector,foldname)
                trainset=other_args.foldlabels(foldsearch).labels(fold).foldname
           end
       end
       end
   else if strcmp(xvalselector,'runs')
      trainset='trainruns';
end
   end
   patheader.vol.fname=[printregressor '_' xvalselector '_' trainset '_srchacc_' num2str(fold)]; %update the name in the header just so things aren't confusing
   newheader=patheader
   newheader.vol={patheader.vol}; %oh matlab
   subj=set_objfield(subj, 'pattern', searchlightmaps{fold}, 'header', newheader); %toss that header into your newly created accuracy pattern 
   
   %duplicate each pattern and set new mat that is diff in accuracy from
   %chance
   subj=duplicate_object(subj, 'pattern', searchlightmaps{fold}, [searchlightmaps{fold} '_minus' num2str(chancelevel)]);
   rawpat=get_mat(subj, 'pattern', searchlightmaps{fold}); %take the existing acc map
   diffpat=rawpat-chancelevel; %scale by chance level for this disc
   subj=set_mat(subj, 'pattern', [searchlightmaps{fold} '_minus' num2str(chancelevel)], diffpat); %set this mat into new pattern objects
   subj=set_objfield(subj,'pattern',[searchlightmaps{fold} '_minus' num2str(chancelevel)], 'group_name',[outputpatname '_minus' num2str(chancelevel)]); %rename group
   %print raw acc and scaled to chance
   save =write_to_analyze(subj,'pattern',searchlightmaps{fold},'output_filename', [printregressor '_' xvalselector '_' trainset '_srchacc_' num2str(fold) other_args.subsetstring], 'pathname', saveimgsdir);
   save =write_to_analyze(subj,'pattern',[searchlightmaps{fold} '_minus' num2str(chancelevel)],'output_filename', [printregressor '_' xvalselector '_' trainset '_srchacc_minus' num2str(chancelevel) '_' num2str(fold) other_args.subsetstring], 'pathname', saveimgsdir);
   patheader.vol
end


end





