function prep_for_mvpaptb_template(study, task, subjectlist, resultfolder, runsrange)
%created by AES 4/16/13, takes data from standard saxelab SPM organization
%and gets it ready for mvpa analysis using princeton toolbox

rootdir=['/mindhive/saxelab2/'];
groupmpvaptbdir='mvpaptb/'
studydir=[rootdir study '/'];
prefix='swrf'; % type of preprocessed image to look for
TR=2;
exwartEvent=0; % 0 to just exclude individual artifact timepoint, 1 to include whole event surrounding artifact
wart='tpoint' % specifies that art exlusion is done on timepoints (rather than events)
responsefunction='hrf';
minVoxels=20; % a given subject's ROI needs at least this many voxels to be used in analyses. if you are going to use rank-based feature selection ultimately, would be most straightforward to just have this number be the same as whatever you decide to use as your other_args.fixednum in run_classification_EIB.m
fextension='.img'; %what kind of images are we looking for?

if ~iscell(subjectlist)
    subjectlist={subjectlist};
end

numSubj=length(subjectlist);
runs=listbolds(task, subjectlist); % calls a simple little script in scripts/aesscripts that finds which bolds correspond to particular task for each subject (by calling a file called EIB_subject_taskruns.mat that lives in the main EIB directory)

for s=1:numSubj

subject=subjectlist{s}

    %figure out what bold dirs this subject has
    boldlist=runs{s}; % get list of bolds for this task for this subject
    boldlist=boldlist(runsrange)    
    runsincluded=[int2str(runsrange(1)) 'to' int2str(runsrange(end))];
    if ~isempty(boldlist) % assuming they have some bold dir for this task, do...

numruns=length(boldlist)    
subjdir=[studydir subject '/'];
cd(subjdir)
mvpadir='mvpa_ptb/';
mkdir(mvpadir)
    
% get whole brain masks (this assumes we are using normalized data)
maskdir='3danat/';
mvpa_mask=dir([maskdir 'skull_strip_mask' fextension]);
mvpa_mask_description{1}='whole_brain_skullstripped';
maskname{1}=[maskdir mvpa_mask.name];

%%%%%%%%%%%%%%%%% MODIFY THESE FOR YOUR HYPOTHESES %%%%%%%%%%%%%%%%%%%%
%hardcode the rois you are interested in

%individual rois stored in each subject's folder
ROISind={'RTPJ_tomloc','LTPJ_tomloc','RSTS_tomloc','LSTS_tomloc','PC_tomloc', 'MMPFC_tomloc', 'DMPFC_tomloc', 'VMPFC_tomloc', 'ROI_rpSTS_BDbiomot', 'ROI_rFFA_kanparcel_EmoBioLoc','ROI_lFFA_kanparcel_EmoBioLoc', 'ROI_rSTS_kanparcel_EmoBioLoc', 'ROI_lSTS_kanparcel_EmoBioLoc', 'ROI_rOFA_kanparcel_EmoBioLoc', 'ROI_lOFA_kanparcel_EmoBioLoc'}; %these refer to actual file names for rois stored in your roi dir
for roin=1:length(ROISind)
    roi_descriptions_ind{roin}=[ROISind{roin} '_ind']; % this is the name by which you refer to the roi (I distinguish between group and ind rois)
end
%group rois that stored elsewhere
ROISgroup={
    %'/mindhive/saxelab/roi_library/functional/EIBrois/ROI_rinsula_wfu_xyz.img',...
    %'/mindhive/saxelab/roi_library/functional/EIBrois/ROI_rvSTR_reward_xyz.img',...
    %'/mindhive/saxelab/roi_library/functional/EIBrois/ROI_vmPFC_reward_xyz.img',...
    %'/mindhive/saxelab/roi_library/functional/EIBrois/ROI_right_ant_temporal_xyz.img',...
    %'/mindhive/saxelab/roi_library/functional/EIBrois/ROI_ramygdala_wfu_xyz.img',...
    '/mindhive/saxelab/roi_library/functional/EIBrois/ROI_MPFC_peelenpeak_xyz.img',...
    '/mindhive/saxelab/roi_library/functional/EIBrois/ROI_lSTS_peelenpeak_xyz.img'
    %'/mindhive/saxelab/roi_library/functional/EIBrois/ROI_rSTS_peelenflip_xyz.img',...
    %'/mindhive/saxelab/roi_library/functional/EIBrois/ROI_lvSTR_reward_xyz.img',...
    %'/mindhive/saxelab/roi_library/functional/EIBrois/ROI_linsula_wfu_xyz.img',...
    %'/mindhive/saxelab/roi_library/functional/EIBrois/ROI_left_ant_temporal_xyz.img',...
    %'/mindhive/saxelab/roi_library/functional/EIBrois/ROI_lamygdala_wfu_xyz.img'
    };

roi_descriptions_group={
    %'rinsula_wfu_xyz_group',...
    %'rvSTR_reward_xyz_group',...
    %'vmPFC_reward_xyz_group',...
    %'right_ant_temporal_xyz_group',...
    %'ramygdala_wfu_xyz_group',...
    'MPFC_peelenpeak_xyz_group',...
    'lSTS_peelenpeak_xyz_group'
    %'rSTS_peelenflip_xyz_group',...
    %'lvSTR_reward_xyz_group',...
    %'linsula_wfu_xyz_group',...
    %'left_ant_temporal_xyz_group',...
    %'lamygdala_wfu_xyz_group'
    };


%roidir=[subjdir 'roi/'];
roidir=[subjdir 'autoROI/'];
numindrois=length(ROISind);
numgrouprois=length(ROISgroup);
count=1;
voxelcount(s,:)=zeros(1,numindrois+numgrouprois); %make a matrix of the number of voxels each subject has in each roi
voxelcountlabel=cell(1,numindrois+numgrouprois); %make an array of names corresponding to each of those rois

% get invidual rois
for roin=1:numindrois
        roiname=ROISind{roin};
        voxelcountlabel{roin}=roiname;
		searchROI=dir([roidir roiname '*' fextension]);
        if ~isempty(searchROI) %if the subject has this roi...
            roi_vol=spm_vol([roidir searchROI.name]);
            roi_mat=spm_read_vols(roi_vol);
            numVoxels=sum(roi_mat(:)); % find out how many voxels are in it
            voxelcount(s, roin)=numVoxels;
            if numVoxels>minVoxels %if it is more than your specified min
            count= count+1;
            mvpa_mask_description{count}=roi_descriptions_ind{roin}; %add it to the list of rois you'll use for this subject
            maskname{count}=[roidir searchROI.name];
            end
        end
end

% now get group rois
for roin=1:numgrouprois
    voxelcountlabel{roin+numindrois}=ROISgroup{roin};
    roi_vol=spm_vol(ROISgroup{roin});
    roi_mat=spm_read_vols(roi_vol);
    numVoxels=sum(roi_mat(:));
    count=count+1;
    maskname{count}=[ROISgroup{roin}];
    mvpa_mask_description{count}=roi_descriptions_group{roin}; 
    voxelcount(s, roin+numindrois)=numVoxels;
end
subjvoxelcount=voxelcount(s,:);


% get some design details from the results folders (this result folder is
% only being used to grab an SPM.mat to figure out the number of timepoints and the onsets, so
% it shouldn't matter what the modeling parameters were. there do need to exist 'art_regression_outliers*.mat')
load(['results/' resultfolder '/SPM.mat'])
numTimepoints=size(SPM.xY.VY,1);
IPS=numTimepoints/numruns; %assumes all runs are same length
ordereddurations=zeros(numTimepoints,2)
%here we are just going to pull up run 1 to get relevant info about number of conditions
%and durations of each 
%(ghetto constraint: the following assumes that the # of conditions in run 1 is the # of conditions in the whole exp)
    r=1;  
    taskStruct=SPM.Sess(r).U;
    numConds=size(taskStruct,2);
%initialize regressors and selectors matrices
    binarizedregressors=zeros(numTimepoints, numConds);
    blocklabelmatrix=zeros(numTimepoints, 1);
    runselectors=zeros(numTimepoints, 1); 
    artselectors=ones(numTimepoints, 1);
    
    for c=1:numConds
        allcond_names{c}=SPM.Sess(r).U(c).name{:}; % these are all the names of conditions in this exp
        numtrialsincond(c)=length(SPM.Sess(r).U(c).dur(1));
    end
    maxtrialsinsinglecond=max(numtrialsincond); %how big dos my trialDurations matrix need to be?
    trialDurations=zeros(c,maxtrialsinsinglecond);
     for c=1:numConds
     trialDurations(c,:)=ceil(SPM.Sess(r).U(c).dur(:)); %% rounds up duration up to nearest TR
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
    runselectors(startRow:endRow,:)=runselectors(startRow:endRow,:)+r; % selector vector specifying the run #
    evenoddselectors(startRow:endRow,:)=evenodd; %creating a selector vector specifying whether we are in an even or odd run

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
        ordereddurations((r-1)*IPS+onset,trialDurations(c,o)); %puts durations in vector ordered by onset (position in the full timecourse)
        ordereddurations=ordereddurations(:,2);
        ordereddurations(ordereddurations==0)=[]; % yields single vector of durations that are ordered by position of trial in timecourse
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
    artselectors(startRow:endRow,1)=artselectors(startRow:endRow,1)-runarts; %this vector is full of ones, so we are going to subtract out the outlier timepoints to make it 0 for every excluded timepoint
    artselectors(artselectors<1)=0; %erm, just in case?
end

% go through  blocklabel vector. if the event has a condition, relabel it
% with a separate number designating the event (this will be used later for
% creating event examples that average across timepoints)
count=0;
blocklabels=blocklabelmatrix;
for x=1:length(blocklabelmatrix)
       value=blocklabelmatrix(x); %this is the condition of this timepoint
       if value~=0 %if the vector has a value...
       count=count+1; %we've come across another example
       duration=ordereddurations(count); %figure out its duration
       blocklabels(x)=count; % h
       blocklabels(x:x+(duration-1))=count; %label that timepoint and subsequent timepoints as a new example
       end
end


%make selector to exclude all rest points (=0 for any rest timepoint)
restselectors=blocklabels;
restselectors(restselectors>0)=1;

%make single selector so exclude all rest AND artifact points
combselectors=zeros(numTimepoints, 1);
combselectors(restselectors==1 & artselectors==1)=1;

%make block selector that also excludes artifact points (already excludes
%rest points)
blocklabelscomb_old=blocklabels;
blocklabelscomb_old(artselectors==0)=0;

if exwartEvent %if excluding whole event
    excluded=blocklabels(artselectors==0); %find event label for all artifacts
    blocklabelscomb_old(blocklabelscomb_old==any(excluded))=0; % sset anytimepoint that is labed as that event to zero
end
blocklabelscomb=blocklabelscomb_old;
%these were labeled based on total number of events, but you've now
%potentially lost some events due to artifacts. relabel as 1 through total number remaining
%events
oldmax=max(blocklabelscomb_old);
count=0;
for c=1:oldmax
    temp=find(blocklabelscomb==c);
    if ~isempty(temp)
       count=count+1;
       blocklabelscomb(blocklabelscomb==c)=count;
    end
end

%%%%%%%%CREATE YOUR SELECTORS%%%%%%%%%%%%%%%%%%%%%%%%%%
% make some selectors for doing cross validation across stimulus dimensions
%check your SPM.mat for the order of your conditions to create a contrast
%like the one below (not a real contrast just a relabeling of conditions)
stimselector= makeselector(binarizedregressors, [1 1 1 1 2 2 2 2]);  %% first fold=stim1 (e.g. auditory), second = stim2 (e.g. visual)
foldlabels(1).name='stimselector';
foldlabels(1).labels(1).foldname='trainstim1';
foldlabels(1).labels(2).foldname='trainstim2';


%% convolve binarized regressor with hrf (umm, currently not using this for anything)
hrfinfo.dt=TR;
hrfinfo.name=responsefunction;
bf = spm_get_bf(hrfinfo);
for c=1:numConds
U.u=binarizedregressors(:,c);
U.name={'reg'}; % just because it needs a name
convolvedregressors(:,c) = spm_Volterra(U, bf.bf);
end


%%%%%%%%MODIFY YOUR REGRESSORS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make regressors for specific descriminations of interest
% first regressor matrix should include single regressor for each
% individual condition
% all remaining regressors should include just binary contrasts (or, if you
% include more than 1 multi cond contrast, make appropriate changes in classification script so
% that you skip all multis. currently assumes that only first is multicond, and skips this when conducting binary classifications)
alleight=reduce_regressors(binarizedregressors, convolvedregressors, [1 2 3 4 5 6 7 8], {'c1','c2','c3','c4','c5','c6','c7','c8'}, 'none', blocklabels, blocklabelscomb); %cell array specifies condition labels, the string e.g. 'none' or 'stimselector' specifies whether there is a non-run partitioning (created above)
collc1VScollc2=reduce_regressors(binarizedregressors, convolvedregressors, [1 1 1 1 2 2 2 2], {'collapsedc1', 'collapsedc2'}, 'stimselector', blocklabels, blocklabelscomb);

%%

%%list of all the bold images to use
boldnames=mat2cell(SPM.xY.P, ones(numTimepoints,1));

%this is just some silliness to ensure intuitive naming that matches with
%other scripts. should just fix these above.
allconds_convolved=convolvedregressors'; %% note these are transposed to be appropriate orientation for the toolbox
allconds_binarized=binarizedregressors'; %% ditto
clearvars convolvedregressors binarizedregressors

%print some stuff just to visualize/sanity check
p= imagesc([blocklabels/max(blocklabels) blocklabelscomb/max(blocklabelscomb) artselectors/max(artselectors) restselectors/max(restselectors) combselectors/max(combselectors) runselectors/max(runselectors) evenoddselectors/max(evenoddselectors) stimselector/max(stimselector) crossmatchedruns/max(crossmatchedruns) contextselector/max(contextselector) genderselector/max(genderselector) crossrunsONEselector/max(crossrunsONEselector) crossrunsTWOselector/max(crossrunsTWOselector)]);
plotselectors=gcf;
saveas(plotselectors, [mvpadir 'selectorvis']);
clear gcf
close all
p=imagesc([(allconds_binarized*.5)' allconds_convolved']);
plotregressors=gcf;
saveas(plotregressors, [mvpadir 'regressorvis']);
clear gcf
close all

%UPDATE WHAT YOU SAVE
%save the relevant structures for the subject
 save([mvpadir 'subjinfo_' task '_' runsincluded], 'boldnames', 'maskname', 'mvpa_mask_description', 'allconds_binarized', 'allconds_convolved',  'allcond_names','blocklabels', 'blocklabelscomb', 'runselectors', 'evenoddselectors', 'artselectors', 'restselectors', 'combselectors', 'stimselector','foldlabels', 'TR', 'subjvoxelcount', 'voxelcountlabel', 'wart');
 save([mvpadir 'discriminations_' task '_' runsincluded], 'alleight', 'collc1VScollc2');
    end
clearvars 'boldnames' 'maskname' 'mvpa_mask_description' 'allconds_binarized' 'allconds_convolved' 'allcond_names' 'blocklabels'  'blocklabelscomb'  'runselectors'  'evenoddselectors' 'artselectors' 'restselectors' 'combselectors' 'stimselector' 'foldlabels'
clearvars 'alleight' 'collc1VScollc2'
end
    voxelcount(voxelcount==0)=nan;
    voxelcount(numSubj+1,:)=nanmean(voxelcount); %add a row with mean # of voxels in the ROI across subjects
    voxelcount(numSubj+2,:)=min(voxelcount); %add a row with min # of voxels in the ROI in any subject
    save([studydir groupmvpaptbdir 'voxelcountsummary.mat'], 'voxelcount', 'voxelcountlabel'); %save a general summary of the rois across subjects
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
