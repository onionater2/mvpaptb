function [subj, results, printregressor, outputthresh]=mvpaptb_classify(rootdir, study, task, runsincluded, mvparootdir, subjectID, subj, condregressors_binarized, condregressors_convolved, condnames, mask_description, xvalselector, averagingselector, averagingselectorcomb, other_args, class_args)
% created by AES 4/18/13
% called by run_classification_EIB
% takes as input a single binary regressor and convolved regressor (need to
% run separate times for separate classifications)
%

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
   printregressor=condregressors_binarized;
else
   shiftedregressor='conds_convolved'; 
   printregressor=condregressors_convolved;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRE-PROCESSING - z-scoring in time and no-peeking anova

% we want to z-score the EPI data (called 'epi'),
% individually on each run (using the 'runs' selectors)
if other_args.zscore
subj = zscore_runs(subj,other_args.bolds,'runs'); %note zscoring is occuring over full timecourse (of run), not just examples
    fslabel='_z';
    if other_args.averaged
    fslabel='_z_avg'; %so that below you grab the right bolds for feature selection
    subj = average_object(subj,'pattern',[other_args.bolds '_z'],averagingselector);
    subj = average_object(subj,'regressors',shiftedregressor,averagingselector);
    subj = average_object(subj,'selector',xvalselector,averagingselector); 
    end
else if ~other_args.zscore
    fslabel='';
    if other_args.averaged
        fslabel='_avg'; %so that below you grab the right bolds for feature selection
    subj = average_object(subj,'pattern',other_args.bolds, averagingselector);
    subj = average_object(subj,'regressors',shiftedregressor,averagingselector);
    subj = average_object(subj,'selector',xvalselector,averagingselector); 
    end
    end
end

% now, create selector indices for the n different iterations of
% the nminusone %ignoring jumbled runs to enable cross validation on just
% even/odd, for example
if other_args.averaged
    crossvalidation_selector=[useselector '_xval'];
    subj = create_xvalid_indices(subj,useselector, 'ignore_jumbled_runs', 'true');
else
    crossvalidation_selector=[useselector '_' additionalselector '_xval'];
    subj = create_xvalid_indices(subj,useselector, 'actives_selname',additionalselector,'new_selstem', crossvalidation_selector, 'ignore_jumbled_runs', 'true');
end

%set things up for feature selection
if strcmp(other_args.fsfunc,'anova')
    fspatterns=[other_args.bolds fslabel];
    fsregressor=finalregressor; %take regressors for contrast of interest
    fsselector=crossvalidation_selector; %take xval selectors for contrast of interest
    selnames = find_group(subj,'selector',fsselector);
else if strcmp(other_args.fsfunc,'activity')
        if other_args.zscore
        fspatterns=[other_args.bolds '_z'];
        else
        fspatterns=other_args.bolds;
        end
        fsregressor=finalregressor;
        fsselector='rests';
        selnames = find_group(subj,'selector',crossvalidation_selector);
    end
end


        % some feature selection functions (e.g. activity) need to know the number of xval folds you'll be using
        % so, find the selectors within the crossvalidiation_selector group
        numfolds = length(selnames);


% run the anova multiple times, separately for each iteration,
% using the selector indices created above
feature_selection_func=str2func(['feature_select_' other_args.fsfunc]);
if ~other_args.featureselect
     numVoxels=sum(subj.masks{1}.mat(:)); % this assumes that first mask is the main ROI 
     xvalmask=mask_description;
else if other_args.featureselect
[subj] = feature_selection_func(subj,fspatterns,fsregressor,fsselector, numfolds, other_args.fixed, other_args.fixednum, 'thresh',other_args.fsthreshold);
 
if ~other_args.fixed
%% this section checks number of voxels surviving at the specified threshold
% iterative process. first take voxels at specified threshold. if insufficient
% number of voxels at that criteria, try again at that threshold + .5
patternTooSmall=1;
while patternTooSmall && other_args.fsthreshold<=1
xvalmask=[fspatterns '_thresh' num2str(other_args.fsthreshold)];
numMasks=length(subj.masks);
numPatterns=length(subj.patterns);
numVoxels=[];
maskIndices=[];
patternIndices=[];

%% find the thresholded masks and count their voxels
for xmask=1:numMasks
    mname=subj.masks{xmask}.name;
    if ~isempty(strfind(mname, 'thresh'))
    maskIndices=[maskIndices xmask];    
    numVoxels=[numVoxels sum(subj.masks{xmask}.mat(:))]; % this assumes that only masks are the main ROI and the xval masks
    end
end

toofew=find(numVoxels<other_args.voxelthreshold);
if isempty(toofew)
    patternTooSmall=0;
else 
    disp(['too few voxels at threshold of ' num2str(other_args.fsthreshold) '. recalculating xval mask with threshold ' num2str(other_args.fsthreshold+.05)])
    other_args.fsthreshold=other_args.fsthreshold+0.05;
    
    % we don't need the masks generated by each iteration of this
    % figure out the masks you want to keep
    allmindices=1:numMasks; 
    for p=1:length(allmindices) %go throug each mask index
        keeperm(p)=~any(allmindices(p)==maskIndices);   
    end
    keepersm=find(keeperm);
    for x=1:length(keepersm)
       temp{x}=subj.masks{x}; 
    end
    subj.masks=temp;
    length(subj.masks);
    temp=[];
    
    % find the previously created feature-selected patterns
    for p=1:numPatterns
       pname=subj.patterns{p}.name;
       if ~isempty(strfind(pname, 'anova'))
           patternIndices=[patternIndices p];
       end
    end
    % figure out the patterns you want to keep
    allpindices=[1:numPatterns];
    for p=1:length(allpindices)
        keeper(p)=~any(allpindices(p)==patternIndices);    
    end
    keepers=find(keeper);
    for x=1:length(keepers)
       temp{x}=subj.patterns{x}; 
    end
    subj.patterns=temp;
    length(subj.patterns);
    temp=[];

    %% rerun the feature selection with more lenient threshold
    [subj] = feature_selection_func(subj,fspatterns,finalregressor,crossvalidation_selector, numfolds, other_args.fixed, other_args.fixednum, 'thresh',other_args.fsthreshold);
     xvalmask=[fspatterns '_thresh' num2str(other_args.fsthreshold)];

end
end
else if other_args.fixed
        xvalmask=[fspatterns '_top' num2str(other_args.fixednum)]; %set this to the group name for the masks
        numMasks=length(subj.masks);
        numVoxels=[];
        for xmask=1:numMasks
        numVoxels=[numVoxels sum(subj.masks{xmask}.mat(:))]; % this assumes that only masks are the main ROI and the xval masks
        end
    end
end
    end
end

regressorForClassifier=finalregressor;
%%% scramble here
if other_args.scramble
    scramname=[finalregressor '_scrambled'];
    subj = scramble_regressors_aes(subj,finalregressor,useselector,scramname, other_args.itersuffix);
    regressorForClassifier=scramname;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLASSIFICATION - n-minus-one cross-validation

% set some basic arguments for the classifier
class_args.train_funct_name = ['train_' other_args.classifier];
class_args.test_funct_name = ['test_' other_args.classifier];
class_args.ignore_1ofn = 'false';
avgnumvoxels=mean(numVoxels);
if exist('class_args.penalty')
    class_args.penalty=class_args.penalty*floor(avgnumvoxels);
end
    
% now, run the classification multiple times, training and testing
% on different subsets of the data on each iteration
disp(['classifying ' subjectID ' for ' printregressor ' in ' mask_description])
[hidethatcraziness, subj, results] = evalc('cross_validation(subj,patterns,regressorForClassifier,crossvalidation_selector,xvalmask,class_args)'); %evalc lets you hide all the normal output of this function and subfunctions

outputthresh=other_args.fsthreshold;
end