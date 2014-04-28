function prep_for_mvpaptb_EIB_main_crosssubj(study, task, inputsubjectlist, resultfolder, runsrange)
%prep_for_mvpaptb_EIB_main('EIB', 'EIB_main', makeIDs('EIB', [28:32]), 'EIB_main_with_art_reg_results_normed', [1:8])
%created by AES 4/16/13, takes data from standard saxelab SPM organization
%and gets it ready for mvpa analysis using princeton toolbox

rootdir=['/mindhive/saxelab2/'];
studydir=[rootdir study '/'];
prefix='swrf'; % type of preprocessed image to look for
TR=2;
responsefunction='hrf';
minVoxels=20; % a given subject's ROI needs at least this many voxels to be used in analyses. if you are going to use rank-based feature selection ultimately, would be most straightforward to just have this number be the same as whatever you decide to use as your other_args.fixednum in run_classification_EIB.m
fextension='.img'; %what kind of images are we looking for?
exwartEvent=0; % 0 to just exclude individual artifact timepoint, 1 to include whole event surrounding artifact
if exwartEvent wart='event'; else wart='tpoint'; end
includeStimLabels=0; % (default=0) will make regressor of specific stimulus labels 
if includeStimLabels
findstims=[rootdir study '/mvpaptb/subject_stim_orders.mat'];%expect stimlabels in a single timecourse in this file (only needed if includeStimLabels=1). will need fields group_itemnumbers, group_itemnames, group_keys (See below)
stims=load(findstims)
end
mvpadir='/mindhive/saxelab2/EIB/crosssubj_mvpa/mvpaptb/'

if ~iscell(inputsubjectlist)
    subjectlist={subjectlist};
end

numSubj=length(inputsubjectlist);
runs=listbolds(task, inputsubjectlist); % calls a simple little script in scripts/aesscripts that finds which bolds correspond to particular task for each subject (by calling a file called EIB_subject_taskruns.mat that lives in the main EIB directory)

% get whole brain masks (this assumes we are using normalized data)
maskdir=['/mindhive/saxelab2/EIB/crosssubj_mvpa/3danat/'];
mvpa_mask=dir([maskdir 'binarized40percent_grey_matter_MNI_fromSPMapriori' fextension]);
mvpa_mask_description{1}='whole_brain_grey_matter';
maskname{1}=[maskdir mvpa_mask.name];


maskname={
    %'/mindhive/saxelab/roi_library/functional/EIBrois/ROI_rinsula_wfu_xyz.img',...
    %'/mindhive/saxelab/roi_library/functional/EIBrois/ROI_rvSTR_reward_xyz.img',...
    %'/mindhive/saxelab/roi_library/functional/EIBrois/ROI_vmPFC_reward_xyz.img',...
    %'/mindhive/saxelab/roi_library/functional/EIBrois/ROI_right_ant_temporal_xyz.img',...
    %'/mindhive/saxelab/roi_library/functional/EIBrois/ROI_ramygdala_wfu_xyz.img',...
    '/mindhive/saxelab/roi_library/functional/EIBrois/ROI_MPFC_peelenpeak_xyz.img',...
    '/mindhive/saxelab/roi_library/functional/EIBrois/ROI_lSTS_peelenpeak_xyz.img',...
    '/mindhive/saxelab/roi_library/functional/EIBrois/ROI_rSTS_peelenflip_xyz.img',...
    '/mindhive/saxelab2/EIB/AllROIsBackup/MPFC_combo.img',...
    '/mindhive/saxelab2/EIB/AllROIsBackup/DMPFC_xyz.img',...
    '/mindhive/saxelab2/EIB/AllROIsBackup/MMPFC_xyz.img',...
    '/mindhive/saxelab2/EIB/AllROIsBackup/VMPFC_xyz.img'
    %'/mindhive/saxelab/roi_library/functional/EIBrois/ROI_lvSTR_reward_xyz.img',...
    %'/mindhive/saxelab/roi_library/functional/EIBrois/ROI_linsula_wfu_xyz.img',...
    %'/mindhive/saxelab/roi_library/functional/EIBrois/ROI_left_ant_temporal_xyz.img',...
    %'/mindhive/saxelab/roi_library/functional/EIBrois/ROI_lamygdala_wfu_xyz.img'
    };

mvpa_mask_description={
    %'rinsula_wfu_xyz_group',...
    %'rvSTR_reward_xyz_group',...
    %'vmPFC_reward_xyz_group',...
    %'right_ant_temporal_xyz_group',...
    %'ramygdala_wfu_xyz_group',...
    'MPFC_peelenpeak_xyz_group',...
    'lSTS_peelenpeak_xyz_group',...
    'rSTS_peelenflip_xyz_group',...
    'MPFC_combo_xyz_group',...
    'DMPFC_tomloc_xyz_group',...
    'MMPFC_tomloc_xyz_group',...
    'VMPFC_tomloc_xyz_group'
    %'lvSTR_reward_xyz_group',...
    %'linsula_wfu_xyz_group',...
    %'left_ant_temporal_xyz_group',...
    %'lamygdala_wfu_xyz_group'
    };


numindrois=0;
numgrouprois=length(maskname);

boldnames=[]
blocklabels=[]  
blocklabelscomb=[]  
runselectors=[]
evenoddselectors=[]
artselectors=[]
restselectors=[]
combselectors=[]
allconds_convolved=[]
allconds_binarized=[]
alleight=[]
faceVScontext=[]
negVSpos=[]
negfVSposf=[]
negcVSposc=[]
maleVSfemale=[]
socialVSnonsoc=[]
socialnVSsocialp=[]
nonsocnVSnonsocp=[]
malenVSmalep=[]
femalenVSfemalep=[]
subjselector=[]
blockcounter=0
blockcombcounter=0
for s=1:numSubj

subject=inputsubjectlist{s}

    %figure out what bold dirs this subject has
    boldlist=runs{s}; % get list of bolds for this task for this subject
    boldlist=boldlist(runsrange)    
    runsincluded=[int2str(runsrange(1)) 'to' int2str(runsrange(end))];
    if ~isempty(boldlist) % assuming they have some bold dir for this task, do...

numruns=length(boldlist)    
subjdir=[studydir subject '/'];
cd(subjdir)

% get some design details from the results folders (this result folder is
% only being used to grab an SPM.mat to figure out the number of timepoints and the onsets, so
% it shouldn't matter what the modeling parameters were. there do need to exist 'art_regression_outliers*.mat')
load(['results/' resultfolder '/SPM.mat'])
numTimepoints=size(SPM.xY.VY,1);
IPS=numTimepoints/numruns; %assumes all runs are same length
ordereddurations=zeros(1,numTimepoints);
%here we are just going to pull up run 1 to get relevant info about number of conditions
%and durations of each 
%(ghetto constraint: the following assumes that the # of conditions in run 1 is the # of conditions in the whole exp, and that the duration of each condition is a constant across events)
    r=1;  
    taskStruct=SPM.Sess(r).U;
    numConds=size(taskStruct,2);
%initialize regressors and selectors matrices
    binarizedregressors=zeros(numTimepoints, numConds);
    blocklabelmatrix=zeros(numTimepoints, 1);
    ind_runselectors=zeros(numTimepoints, 1); 
    ind_artselectors=ones(numTimepoints, 1);
    
    for c=1:numConds
        allcond_names{c}=SPM.Sess(r).U(c).name{:}; % these are all the names of conditions in this exp
        numtrialsincond(c)=length(SPM.Sess(r).U(c).dur); %% assumes each event for a condition is same length across presentations, and rounds up duration up to nearest TR
    end
    maxtrialsinsinglecond=max(numtrialsincond);
    trialDurations=zeros(c,maxtrialsinsinglecond);
     for c=1:numConds
     trialDurations(c,:)=ceil(SPM.Sess(r).U(c).dur(:)); %% code below assumes each event for a condition is same length across presentations, and rounds up duration up to nearest TR
     end

% now go through each run for this task    
for r=1:numruns
    run=boldlist(r);
    if mod(r,2)==0
        evenodd=2; %even
    else
        evenodd=1; %odd
    end
    %we are making one long design matrix but creating it run by run so we
    %need to know where in the matrix we should be adding values
    startRow=(r-1)*numTimepoints/numruns+1;
    endRow=(r-1)*numTimepoints/numruns+numTimepoints/numruns;
    binarizedrunmatrix=zeros(IPS,numConds); %matrix for a single run
    blockrunmatrix=zeros(IPS,1); % single vector that marks each event in the run with number corresponding to condition
    ind_runselectors(startRow:endRow,:)=ind_runselectors(startRow:endRow,:)+r+8*(s-1); % selector vector specifying the run #
    ind_evenoddselectors(startRow:endRow,:)=evenodd; %creating a selector vector specifying whether we are in an even or odd run

% go through each condition and find its duration (assumed to be constant) and its onsets    
for c=1:numConds
    condonsets=SPM.Sess(r).U(c).ons; %rounding assumes that all onsets are on TR. change this f ever jittering off TR
    numonsets=length(condonsets);
    binarizedrunmatrix(condonsets,c)=1; %for each onset, add a 1 to the appropriate condition column
    blockrunmatrix(condonsets,1)=c; %for each onset specify the condition
    
    %right now we only have a 1 at the onset of each event. go through and
    %add 1s to subsequent timepoints based on duration value
    for o=1:numonsets
        onset=condonsets(o);
        binarizedrunmatrix(onset:onset+(trialDurations(c,o)-1),c)=1;
        ordereddurations((r-1)*IPS+onset)=trialDurations(c,o); %puts durations in vector ordered by onset (position in the full timecourse)
    end
end

binarizedregressors(startRow:endRow,:)=binarizedrunmatrix; %plop that single run matrix into the larger design matrix
blocklabelmatrix(startRow:endRow,1)=blockrunmatrix; %plot the single run vector into the whole timecourse vector

    % load art regressors and use to create art_selector vector (=0 for any
    % artifact timepoint)
    if ~isempty(dir(['bold/0' num2str(run) '/art_regression_outliers_sw*.mat']))
    matfile=adir(['bold/0' num2str(run) '/art_regression_outliers_sw*.mat']);
    load(matfile{1})
    runarts=R(:,:); %assumes that all columns are single timepoint exclusion vectors
    else if ~isempty(dir(['bold/0' num2str(run) '/art_regression_outliers_and_movement*sw*.mat']))
    matfile=adir(['bold/0' num2str(run) '/art_regression_outliers_and_movement*sw*.mat']);
    load(matfile{1})
    runarts=R(:,1:end-6); %assumes that the last six columns in this file are motion parameters and that the rest are single timepoint exclusion vectors
    else
       error('it looks like you do not have an art_regression_outliers file in your bold directory');
    end
    end
    runarts=sum(runarts,2); %sum across the regressors to get one vector with a 1 for every outlier timepoint
    ind_artselectors(startRow:endRow,1)=ind_artselectors(startRow:endRow,1)-runarts; %this vector is full of ones, so we are going to subtract out the outlier timepoints to make it 0 for every excluded timepoint
    ind_artselectors(ind_artselectors<1)=0; %erm, just in case?
end

% go through  blocklabel vector. if the event has a condition, relabel it
% with a separate number designating the event (this will be used later for
% creating event examples that average across timepoints)
count=blockcounter;
localcount=0
ind_blocklabels=blocklabelmatrix;
ordereddurations(ordereddurations==0)=[]; % yields single vector of durations that are ordered by position of trial in timecourse   
for x=1:length(blocklabelmatrix)
       value=blocklabelmatrix(x); %this is the condition of this timepoint
       if value~=0 %if the vector has a value...
       count=count+1; %we've come across another example
       localcount=localcount+1
       duration=ordereddurations(localcount); %figure out its duration
       ind_blocklabels(x)=count; % h
       ind_blocklabels(x:x+(duration-1))=count; %label that timepoint and subsequent timepoints as a new example
       end
end
blockcounter=max(ind_blocklabels) %update so that you begin counting with next value for next subject


%make selector to exclude all rest points (=0 for any rest timepoint)
ind_restselectors=ind_blocklabels;
ind_restselectors(ind_restselectors>0)=1;

%make single selector so exclude all rest AND artifact points
ind_combselectors=zeros(numTimepoints, 1);
ind_combselectors(ind_restselectors==1 & ind_artselectors==1)=1;

%make block selector that also excludes artifact points (already excludes
%rest points)
blocklabelscomb_old=ind_blocklabels;
blocklabelscomb_old(ind_artselectors==0)=0;
if exwartEvent %if excluding whole event
    excluded=blocklabels(artselectors==0); %find event label for all artifacts
    blocklabelscomb_old(blocklabelscomb_old==any(excluded))=0; % sset anytimepoint that is labed as that event to zero
end
ind_blocklabelscomb=blocklabelscomb_old;
%these were labeled based on total number of events, but you've now
%potentially lost some events due to artifacts. relabel as 1 through total number remaining
%events
oldmax=max(blocklabelscomb_old);
count=blockcombcounter;
for c=1:oldmax
    temp=find(ind_blocklabelscomb==c);
    if ~isempty(temp)
       count=count+1;
       ind_blocklabelscomb(ind_blocklabelscomb==c)=count;
    end
end
blockcombcounter=max(ind_blocklabelscomb) %update so that you begin counting with next value for next subject


%% convolve binarized regressor with hrf (umm, currently not using this for anything)
hrfinfo.dt=TR;
hrfinfo.name=responsefunction;
bf = spm_get_bf(hrfinfo);
for c=1:numConds
U.u=binarizedregressors(:,c);
U.name={'reg'}; % just because it needs a name
convolvedregressors(:,c) = spm_Volterra(U, bf.bf);
end

%%

%%list of all the bold images to use
ind_boldnames=mat2cell(SPM.xY.P, ones(numTimepoints,1)); %%%append all subjs

%this is just some silliness to ensure intuitive naming that matches with
%other scripts. should just fix these above.
allconds_convolved=[allconds_convolved;convolvedregressors]; %% note these are transposed to be appropriate orientation for the toolbox
allconds_binarized=[allconds_binarized;binarizedregressors]; %% ditto
clearvars convolvedregressors binarizedregressors
blocklabels=[blocklabels; ind_blocklabels]  
blocklabelscomb=[blocklabelscomb;ind_blocklabelscomb]  
runselectors=[runselectors;ind_runselectors]
evenoddselectors=[evenoddselectors;ind_evenoddselectors]
artselectors=[artselectors;ind_artselectors]
restselectors=[restselectors;ind_restselectors]
combselectors=[combselectors;ind_combselectors]
subjselector=[subjselector;ind_runselectors*0+s]
boldnames=[boldnames;ind_boldnames]
    
clearvars 'ind_boldnames' 'ind_blocklabels'  'ind_blocklabelscomb'  'ind_runselectors'  'ind_evenoddselectors' 'ind_artselectors' 'ind_restselectors' 'ind_combselectors'
clearvars 'ind_alleight' 'ind_faceVScontext' 'ind_negVSpos' 'ind_negfVSposf' 'ind_negcVSposc' 'ind_maleVSfemale' 'ind_socialVSnonsoc' 'ind_socialnVSsocialp' 'ind_nonsocnVSnonsocp' 'ind_malenVSmalep' 'ind_femalenVSfemalep'
    
    end
end

%make regressors for specific descriminations of interest
%by using blocklabelscomb we are ignoring any example that has a single
%timepoint missing as an art outlier
% first regressor matrix should include single regressor for each
% individual condition
% all remaining regressors should include just binary contrasts (or, if
% more than 1 multi cond contrast, make appropriate changes in EIB_main so
% that you skip all multis)
alleight=reduce_regressors(allconds_binarized, allconds_convolved, [1 2 3 4 5 6 7 8], {'mu','fu','mh','fh','nu','su','nh','sh'}, 'subjselector', blocklabels, blocklabelscomb);
faceVScontext=reduce_regressors(allconds_binarized, allconds_convolved, [1 1 1 1 2 2 2 2], {'face', 'context'}, 'subjselector', blocklabels, blocklabelscomb);
negVSpos=reduce_regressors(allconds_binarized, allconds_convolved, [1 1 2 2 1 1 2 2], {'neg', 'pos'}, 'subjselector', blocklabels, blocklabelscomb);
negfVSposf=reduce_regressors(allconds_binarized, allconds_convolved, [1 1 2 2 0 0 0 0], {'negf', 'posf'}, 'subjselector', blocklabels, blocklabelscomb);
negcVSposc=reduce_regressors(allconds_binarized, allconds_convolved, [0 0 0 0 1 1 2 2], {'negc', 'posc'}, 'subjselector', blocklabels, blocklabelscomb);
maleVSfemale=reduce_regressors(allconds_binarized, allconds_convolved, [1 2 1 2 0 0 0 0], {'male', 'female'}, 'subjselector', blocklabels, blocklabelscomb);
socialVSnonsoc=reduce_regressors(allconds_binarized, allconds_convolved, [0 0 0 0 1 2 1 2], {'nonsoc', 'social'}, 'subjselector', blocklabels, blocklabelscomb);
socialnVSsocialp=reduce_regressors(allconds_binarized, allconds_convolved, [0 0 0 0 0 1 0 2], {'socialn', 'socialp'}, 'subjselector', blocklabels, blocklabelscomb);
nonsocnVSnonsocp=reduce_regressors(allconds_binarized, allconds_convolved, [0 0 0 0 1 0 2 0], {'nonsocn', 'nonsocp'}, 'subjselector', blocklabels, blocklabelscomb);
malenVSmalep=reduce_regressors(allconds_binarized, allconds_convolved, [1 0 2 0 0 0 0 0], {'malen', 'malep'}, 'subjselector', blocklabels, blocklabelscomb);
femalenVSfemalep=reduce_regressors(allconds_binarized, allconds_convolved, [0 1 0 2 0 0 0 0], {'femalen', 'femalep'}, 'subjselector', blocklabels, blocklabelscomb);


allconds_convolved=allconds_convolved'
allconds_binarized=allconds_binarized'
%print some stuff just to visualize/sanity check
p= imagesc([blocklabels/max(blocklabels) blocklabelscomb/max(blocklabelscomb) artselectors/max(artselectors) restselectors/max(restselectors) combselectors/max(combselectors) runselectors/max(runselectors) evenoddselectors/max(evenoddselectors) subjselector/max(subjselector)]);
plotselectors=gcf;
saveas(plotselectors, [mvpadir 'selectorvis']);
clear gcf
close all
p=imagesc([(allconds_binarized*.5)' allconds_convolved']);
plotregressors=gcf;
saveas(plotregressors, [mvpadir 'regressorvis']);
clear gcf
close all

%save the relevant structures for the subject
save([mvpadir 'subjinfo_' task '_' runsincluded], 'boldnames', 'maskname', 'mvpa_mask_description', 'allconds_binarized', 'allconds_convolved',  'allcond_names','blocklabels', 'blocklabelscomb', 'runselectors', 'evenoddselectors', 'artselectors', 'restselectors', 'combselectors', 'subjselector','TR', 'wart');
 save([mvpadir 'discriminations_' task '_' runsincluded], 'alleight', 'faceVScontext', 'negVSpos', 'negfVSposf', 'negcVSposc', 'maleVSfemale', 'socialVSnonsoc','socialnVSsocialp', 'nonsocnVSnonsocp', 'malenVSmalep', 'femalenVSfemalep');
 
end

function [reduced]=reduce_regressors(bregressors, cregressors, newcondindices, newcondnames, acrosssel,blocklabel, blocklabelscomb)
%this takes the main regressor and makes a regressor expressing the
%descrimination of interest
reducedNum=length(newcondnames);
reduced.names=newcondnames;

for cond=1:reducedNum
   indices=find(newcondindices==cond); %figure out all the columns that correspond to conditions that will be part of your "reduced" condition
   thiscondregressors=bregressors(:,indices); 
   thiscondconvolved=cregressors(:,indices);
   summedregressors=sum(thiscondregressors,2); %sum across those columns to get one vector for the new condition
   summedconvolved=sum(thiscondconvolved,2);
   reg(:,cond)=summedregressors;
   regc(:,cond)=summedconvolved;
end

isanexample=sum(reg,2); % 1 for every timepoint that is part of an example in this discrimination
reducedlabel=blocklabel.*isanexample;
oldmax=max(blocklabel); %how many examples were in the full dataset
count=0;
%go through every old example and see if it is in the current set of
%examples
for c=1:oldmax
    temp=find(reducedlabel==c);
    if ~isempty(temp) %when you find one, relabel it with a new example counter
       count=count+1;
       reducedlabel(reducedlabel==c)=count;
    end
end

%do the same thing but for the combblock selector
reducedcomblabel=blocklabelscomb.*isanexample;
oldmax=max(blocklabelscomb);
count=0;
for c=1:oldmax
    temp=find(reducedcomblabel==c);
    if ~isempty(temp)
       count=count+1;
       reducedcomblabel(reducedcomblabel==c)=count;
    end
end


reduced.binarizedreg=reg';
reduced.convolvedreg=regc';
reduced.crossselector=acrosssel;
reduced.averaginglabels=reducedlabel;
reduced.averaginglabelscomb=reducedcomblabel;
end

function selector= makeselector(bregressors, labels)
    numlabels=max(labels);
    selector=0*sum(bregressors,2); %make vector of zeros of same length as regressors
    for lb=1:numlabels
       indices=find(labels==lb); % get indices of all conditions corresponding to this label
       thisselector=bregressors(:,indices);
       summedselector=sum(thisselector,2); %sum together regressors corresponding to those conditions
       selector(summedselector==1)=lb; %label selector vector appropriately
    end
end
