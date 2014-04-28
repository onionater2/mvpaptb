function run_classification_template(study, task, subjectlist, runsincluded, varargin)
%e.g. run_classification_template('EIB', 'EIB_main', makeIDs('EIB', [1:5, 7:9], '1to8', {'folderprefix', 'version2'})
% created by AES 4/16/2013 based on princton mvpa toolbox documentation (see tutorial)
% run startstuff.m from EIB to set paths for this and related analyses
%this executes basic classification of conditions in ROIs
% assumes you have already created a subjinfo.mat file for each subject,
% stored in subject's mvpa_ptb subdirectory (this is created using
% prep_for_mvpaptb)
% main function called is mvpaptb_classify, also found in
% scripts/aesscripts/mvpaptb
% this is intended as an EIB specific wrapper, but shouldn't be hard to
% extend

rootdir=['/mindhive/saxelab2/'];
studydir=[rootdir study '/'];
cd(studydir)
startstuff

%parameters that require no thought)
folderprefix='';
masksubset={};
discsubset={};
studydir=[rootdir study '/'];
printandsave_script='/mindhive/saxelab/scripts/aesscripts/mvpaptb/printandsave_mvpaptb';
uselotsofspace=0; %set this to 1 if you want to be a space gobbler (and get some sanity check figs). save this for when debugging a new analysis/feature
checkbeforedelete=0;
mvparootdir='mvpa_ptb/'; %specify folder containing individual subject info... contains subjinfo with  discriminations and correct brain mask (all analyses prior to 6/7/13 were accidently making whole brain mask from whole brain structural rather than gray matter mask. wooops, should rerun any analyses that use global normalization using graymatter mask?)
other_args.bolds='epi';
other_args.imagetype='swrf'; % check what kind of boldnames you have in subjinfo
%parameters that require some thought
perm_its=1; %if doing permutation analysis, will generate this many iterations (default==1)
scrambled_regs=0; % 1 if scrambled, 0 if real (Default=0)
suffix=''; %if doing permutation analysis, suffix will be added to classification .mat files specifying the iteration

%optional input arguments in cell array e.g. {'rootdir', '/mindhive/saxelab/',
%'scramble', 1}
if ~isempty(varargin)
for vcheck=1:2:length(varargin{1})
   if strcmp(varargin{1}(vcheck), 'rootdir')
       rootdir=varargin{1}{vcheck+1} %number of iterations for permutation analysis
   end 
   if strcmp(varargin{1}(vcheck), 'perm_its')
       perm_its=varargin{1}{vcheck+1} %number of iterations for permutation analysis
   end
   if strcmp(varargin{1}(vcheck), 'scramble')
       scrambled_regs=varargin{1}{vcheck+1}; % 1 if scrambled, 0 if real
   end
   if strcmp(varargin{1}(vcheck), 'folderprefix')
       folderprefix=varargin{1}{vcheck+1};
   end
   if strcmp(varargin{1}(vcheck), 'discsubset')
       discsubset=varargin{1}{vcheck+1}; % cell array containing discriminations of interest
   end
   if strcmp(varargin{1}(vcheck), 'roisubset')
       masksubset=varargin{1}{vcheck+1}; %cell array containing rois of interest
   end
   if strcmp(varargin{1}(vcheck), 'suffix')
       suffix=varargin{1}{vcheck+1}; %suffix that will be added to 
   end
end
end
other_args.scramble=scrambled_regs;
other_args.perm_iterations=perm_its;
other_args.itersuffix=suffix;
other_args.hemodynamic_delay=4; %6 would also be reasonable...
startmask=2;% 2 to skip whole brain mask, 1 to include whole brain (maybe useful if you want to project back informative voxels)
other_args.featureselect=1; %1= feature select using anova, 0= use all voxels in mask
other_args.fixed=1; %alternative=0. this specifies whether you are using rank (fixed) or threshold based feature selection
other_args.fixednum=80; %matters only if fixed=1this specifies the number of voxels to take if using fixed rank-based feature selection
other_args.fsthreshold=0.05; % matters only if fixed=0... this specifies (initial) p-value for thresholding stat maps from anova for each xval fold
other_args.voxelthreshold=80; %matters only if fixed=0... freak out if xval fold has fewer than this many voxels 
other_args.fsfunc='taskvsrest'; %feature selection metric. default='taskvsrest', alternative:'discriminability'
other_args.hpfilter=1; %default=0
other_args.detrend=1; %default=0
other_args.globalnorm=0; %%default=0 (subtract from each timepoint global mean across voxels.)
other_args.zscore=1; %default=1;
%%%  pattern, regressors and selectors to use in analysis
other_args.binary=1;% 0 if using convolved. currently not an option.
other_args.averaged=1; % 0 if using each timepoint. don't have another choice ATM
other_args.wart=1; % 0 if not using artifact regressors
other_args.notes='blah lbah'; % write yourself some  notes about this analysis
%parameters that require more thought)
other_args.classifier='libsvm' % 'lda' 'ridge', 'gnb'
class=other_args.classifier; %for naming folder. simplifed from original naming convention (e.g. bpbp0)
multi=0; % 1 if multiclass (nn w/ back prop, smlr), 0 if only supports binary classification (SVM, logreg)
if strcmp(other_args.classifier, 'gnb') % probably don't want to use unless you've done some preprocessing/dimensionality reduction to ensure independent features. though, it's fast and reaches assymptote faster that other methods (maybe good if you have few examples: see http://ai.stanford.edu/~ang/papers/nips01-discriminativegenerative.pdf)
    class_args.search4c=2; % just specifies that no cost parameter is used here (optimized or fixed)
    multi=0;
else if strcmp(other_args.classifier, 'bp') %see train_bp.m in mvpa/core/learn
    class_args.nHidden =0; %specifies number of hidden layers
    class_args.search4c=2; % just specifies that no cost parameter is used here (optimized or fixed)
    multi=1;
else if strcmp(other_args.classifier, 'logreg') || strcmp(other_args.classifier, 'ridge') %see train_logreg.m in mvpa/core/learn
    class_args.penalty = 0.05; % will multiply this constant times the number of voxels in the last xval fold (the more voxels the more you want to penalize).     %penalty adopted from EBC tutorial: "to select a penalty parameter, we choose the value taken from default_params, which in this case is 0.05, multiplied by the number of voxels included int the analysis (400). A good rule of thumb for ridge regression is that the penalty parameter should increase linearly as the number of voxels included increases."   
    class_args.search4c=0; %no searching, see above
    multi=0;
else if strcmp(other_args.classifier, 'smlr') % sparse logistic regression
    class_args.tol=exp(-3); % this is the default optimization tolerance in smlr script, just including here for clarity
    class_args.lambda=1; % this indicates amount of regularization (default in smlr, but option to change here)
    %see smlr.m to implement greater degree of parameter flexibility here. or blindly accept defaults :)
    class_args.search4c=0;
    multi=1;
    else if strcmp(other_args.classifier, 'libsvm') %see train_libsvm.m in aesscripts/mvpaptb... uses libsvm library stored in aesscripts/svm, downloaded from here: http://www.csie.ntu.edu.tw/~cjlin/libsvm/
    class_args.kernel_type = 0; %0 -- linear: u'*v; 1 -- polynomial: (gamma*u'*v + coef0)^degree; 2 -- radial basis function: exp(-gamma*|u-v|^2); 3 -- sigmoid: tanh(gamma*u'*v + coef0)  
    class_args.cost = 1; %if not searching for optimal c, set a fixed default. people seem happy with 1.
    class_args.search_for_c = 0; % 1= search for optimal c by xvalidating within training set, 0 = uses .cost
        if class_args.search_for_c
            class_args.search4c=1; %  specifies that  cost parameter is optimized within the training set. awkward redundacny...
            class_args.k_fold_xval      = 8; %how many xval folds to do within training set to find optimal C parameter
        else
            class_args.search4c=0; %cost parameter is fixed
        end
    multi=0;
    class_args.svm_type = 0; %svm_type : set type of SVM (default 0): 0 -- C-SVC; 1 -- nu-SVC; 2 -- one-class SVM; 3 -- epsilon-SVR; 4 -- nu-SVR
else if strcmp(other_args.classifier, 'liblinear') 
    multi=0;
    class_args.search4c=0;
    class_args.cost = 1; %unclear why liblinear doesn't have built in search_for_c option like libsvm does? should be easy to extend but need to look into it
    class_args.svm_type = 3; %3; %     set type of solver (default: svm_type=3: L1-loss and L2-regularization) for classification;
    %see train_liblinear for AES notes about this, and key for specifying
    %desired classifier
        end
        end
        end
    end
    end
end
if strcmp(other_args.fsfunc,'taskvsrest')
    other_args.fsfunc='activity'
else if strcmp(other_args.fsfunc,'discriminability')
    other_args.fsfunc='anova'
    end
end

% TO DO: implement multiclass SVM. decide on one vs. all or all one vs. ones.

numSubj=length(subjectlist);
if ~iscell(subjectlist)
    subjectlist={subjectlist};
end

for s=1:numSubj
    subjectID=subjectlist{s}
    
subjdir=[studydir subjectID '/'];
cd(subjdir)
load([mvparootdir 'subjinfo_' task '_' runsincluded]);

disc=load([mvparootdir 'discriminations_' task '_' runsincluded]);
load([mvparootdir 'discriminations_' task '_' runsincluded]);
discnames=fieldnames(disc)
numDisc=length(discnames);
writtenParamFile=0

%if you want to global normalize your patterns later, calculate the global mean timecourse now
if other_args.globalnorm
    % start by creating an empty subj structure
    subj = init_subj(study, subjectID);
    mask=maskname{1}; % whole brain mask
    mask_description=mvpa_mask_description{1};
    %%% create the mask that will be used when loading in the data
    subj = load_analyze_mask(subj,mask_description,mask);
    % now, read and set up the actual data, keeping only the voxels active in the
    % mask (see above)
    temp_args.globalnorm=0; %so that it doesn't try to global normalize here when the global timecourse does not yet exist...
    subj = load_analyze_pattern_AES(subj,other_args.bolds,mask_description, boldnames, temp_args);
    patternsize=subj.patterns{1}.matsize;
    tclength=patternsize(2);
    numvox=patternsize(1);
    for t=1:tclength
        globalmean(t)=mean(subj.patterns{1}.mat(:,t));
    end
    other_args.globaltimecourse=globalmean;
end

%%%%%%%%%%%%
% some details for naming folders.
if other_args.binary boc='bin'; else boc='conv'; end
if other_args.averaged aon='avgd'; else aon='tpts'; end
if other_args.wart won='wart'; else won='noart'; end
if other_args.featureselect fs=['fs' other_args.fsfunc]; else fs='wholemask'; end
if other_args.hpfilter hp='hpfilt'; else hp='nofilt'; end 
if other_args.detrend dt='detr'; else dt='nodet'; end 
if other_args.globalnorm gn='glnormed'; else gn='noglnorm'; end
if other_args.zscore zs='zsc'; else zs='noz'; end
if other_args.fixed fsf=['ranktop' num2str(other_args.fixednum)]; else fsf=['thresh' num2str(other_args.fsthreshold)]; end
hrfshift=['hrfshift' num2str(other_args.hemodynamic_delay)];
if class_args.search4c==1 
    s4c='costopt'; 
else if class_args.search4c==0
        s4c='costfix';
        if exist('class_args.cost')
        s4c=['costfix' class_args.cost]; 
        else if exist('class_args.penalty')
        s4c=['costfixedvoxX' class_args.penalty]; 
            end
        end
else s4c='costspecNA'; 
    end
end

foldername=[folderprefix 'crossruns_newROIS_' task '_' runsincluded '_' class '_' other_args.imagetype '_' boc '_' won '_' hrfshift '_' fs '_' fsf '_' aon '_' zs '_' hp '_' dt '_' gn '_' s4c]
if other_args.scramble
foldername=['null_' foldername]
end
if exist([mvparootdir foldername])
if checkbeforedelete    
removeDir = questdlg('existing results directory found! overwrite?','Dir with same name found','Yes','No','No');
if strcmpi(removeDir,'No');
   error('directory already exists and was not overwritten. rename your old analysis')
else 
    display(['overwriting ' foldername])
end
end
end

numMasks=length(maskname);
minmask=startmask;  % 2 to skip whole brain mask, 1 to include whole brain (maybe useful if you want to project back informative voxels)
maxmask=numMasks;
for m=minmask:maxmask 
%for m=1
    percentcomplete=m/numMasks*100; %just so that you can know how far along you are in a subject when being impatient
    disp(['classified ' num2str(percentcomplete) ' % of ROIs'])


% start by creating an empty subj structure (note: we'll mostly be using these mvpaptb
% functions (init.., get...) for accessing the structure as an object even though it's just a regular structure
% with normal indexing capabilities. certain ghetto hacks violate this convention, however, out of laziness/disregard for beauty)

subj = init_subj(study, subjectID);

mask=maskname{m};
mask_description=mvpa_mask_description{m};
mvpadir = [mvparootdir foldername '/'];
mkdir(mvpadir)
savemasks{m}=[]
subjectsave.masks{m}=[];

%enables running this on only a subset of masks
dothismask=1;
if ~isempty(masksubset)
if any(strcmp(mask_description, masksubset))
    dothismask=1;
else
    dothismask=0;
end 
end

if dothismask    

%%% create the mask that will be used when loading in the data
subj = load_analyze_mask(subj,mask_description,mask);

% now, read and set up the actual data, keeping only the voxels active in the
% mask (see above)
subj = load_analyze_pattern_AES(subj,other_args.bolds,mask_description, boldnames, other_args);
nummaskstemp=size(subj.masks,1);
for mi=1:nummaskstemp
    if strcmp(subj.masks{mi}.name, mask_description)
    savemasks{m}=subj.masks{mi};   
end
end


if multi
    startd=1; %if multiclass you can use all discriminations
else
    startd=2; % else skip the first and use binary only (NOTE: this logic assumes that only first discrimination (all conditions) is multiclass. make this better.)
end
for d=startd:length(discnames)
%for d=6 %temp just want posVSneg
    n=discnames{d}

    %enables running this on only a subset of masks
dothisdisc=1;
if ~isempty(discsubset)
if any(strcmp(n, discsubset))
    dothisdisc=1;
else
    dothisdisc=0;
end 
end
    
if dothisdisc    
xvalselector='runs'; %%assumes that all classifications can be done across runs. some will additionally be done with an across stimulus xval
binarized=[n '.binarizedreg'];
convolved=[n '.convolvedreg']; %currently don't use convolved for anything, but there it is..
condnames=[n '.names'];
across=eval([n '.crossselector']);
averagingselector=[n '.averaginglabels'];
averagingselectorcomb=[n '.averaginglabelscomb'];


if writtenParamFile==0
% print classifier parameters to a text file for when you inevitably forget everything
f=fopen([subjdir mvpadir 'classparams.txt'],'w');

names=fieldnames(other_args);
numnames=size(names, 1);
for x=1:numnames;
    variablename=['other_args.' names{x}];
    variable=eval(variablename);
    if ~ischar(variable)
        variable=num2str(variable);
    end
    fprintf(f, '%s %s', [names{x} ': '], variable);
    fprintf(f,'\n');
end

names=fieldnames(class_args);
numnames=size(names, 1);
for x=1:numnames;
    variablename=['class_args.' names{x}];
    variable=eval(variablename);
    if ~ischar(variable)
        variable=num2str(variable);
    end
    fprintf(f, '%s %s', [names{x} ': '], variable);
    fprintf(f,'\n');
end

    fclose(f);
    
writtenParamFile=1;    
end
%putting this in after printing classparam.txt    
other_args.allconds_binarized=eval('allconds_binarized'); 
    

%classify across runs
disp(['xval selector: ' xvalselector])
disp(['mask (' num2str(m/numMasks*100) '%): ' mask_description])
disp(['discrimination: ' n])

for permrep=1:other_args.perm_iterations
xvalselector='runs'; %%assumes that all classifications can be done across runs. some will additionally be done with an across stimulus xval
[subjectsave, results, printregressor, fsoutputthreshold]=mvpaptb_classify(rootdir, study, task, runsincluded, mvparootdir, subjectID, subj, binarized, convolved, condnames, mask_description, xvalselector, averagingselector, averagingselectorcomb, other_args, class_args);

%want to save everything but the patterns, because those are huge and nick will hate you. go in
%and kill the actual .mats for each pattern
numpatterns=size(subj.patterns);
for p=1:numpatterns
   subjectsave.patterns{p}.mat=0; 
end
directory=pwd;
    if other_args.scramble 
        if strcmp(other_args.itersuffix, '')
    save([mvpadir mask_description '_' printregressor '_' xvalselector '_classification_iter' num2str(permrep) '.mat'], 'results');
        else
    save([mvpadir mask_description '_' printregressor '_' xvalselector '_classification_iter' num2str(other_args.itersuffix) '.mat'], 'results');
        end  
    else
    run(printandsave_script);
    end
    
%% now classify again across other xval folds, if specified (this really should just call a single function with new parameters rather than running through copied + pasted code. fix this asap)
if ~strcmp(across,'none')
for other=1:length(across)%go through all cross selectors 
    if iscell(across)
    xvalselector=across{other};
    else
    xvalselector=across;
    end
    disp(['discrimination: ' n])
    disp(['xval selector: ' xvalselector])
    disp(['mask (' num2str(m/numMasks*100) '%): ' mask_description])
    [subjectsave, results, printregressor, fsoutputthreshold]=mvpaptb_classify(rootdir, study, task, runsincluded, mvparootdir, subjectID, subj, binarized, convolved, condnames, mask_description, xvalselector, averagingselector, averagingselectorcomb, other_args, class_args);
    
    %want to save everything but the patterns, because those are huge. go in
    %and kill the actual .mats for each pattern
    numpatterns=size(subj.patterns);
    for p=1:numpatterns
        subjectsave.patterns{p}.mat=0; 
    end
    directory=pwd;
    if other_args.scramble 
        if strcmp(other_args.itersuffix, '')
    save([mvpadir mask_description '_' printregressor '_' xvalselector '_classification_iter' num2str(permrep) '.mat'], 'results');
        else
    save([mvpadir mask_description '_' printregressor '_' xvalselector '_classification_iter' num2str(other_args.itersuffix) '.mat'], 'results');
        end  
    else
    run(printandsave_script);
    end
    
    
    
    
end   
end
end
end
end
end
end
subj=0; %done with that subject. kill it and move on.
end
savemasks{1}=subjectsave.masks{1};
save([mvpadir  'allmasks.mat'], 'savemasks');
disp('we finished!')
end





function [subj] = load_analyze_pattern_AES(subj,new_patname,maskname,filenames, other_args, varargin)
% minor AES tweaks
%1)keep from listing every volume loaded in the command
%line, because that's just annoying
% 2) added option to global normalize (therefore now takes other_args as an argument)

% Loads an ANALYZE dataset into a subject structure
%
% [SUBJ] = LOAD_ANALYZE_PATTERN(SUBJ,NEW_PATNAME,MASKNAME,FILENAMES,...)
%
% Adds the following objects:
% - pattern object called NEW_PATNAME masked by MASKNAME
%% Options
% NEW_PATNAME is the name of the pattern to be created
% 
% MASKNAME is an existing boolean mask in the same reference space
% that filters which voxels get loaded in. It should 
%
% All patterns need a 'masked_by' mask to be associated with. The mask
% contains information about where the voxels are in the brain, and
% allows two patterns with different subsets of voxels from the same
% reference space to be compared
%
% See the Howtos (xxx) section for tips on loading data without a
% mask
%
% FILENAMES is a cell array of strings, of .IMG filenames to load
% in. Just the stem, not the extension. If FILENAMES is a string,
% it will automatically get turned into a single-cell array for
% you. If the string contains an asterisk, the string will be
% converted into a cell array of all matching files.
%
% e.g. to load in mydata.IMG:
%   subj = load_analyze_pattern(subj,'epi','wholebrain',{'mydata.img'});
%
% SINGLE (optional, default = false). If true, will store the data
% as singles, rather than doubles, to save memory. Until recently,
% Matlab could store as singles, but none of its functions could do
% much with them. That's been improved, but it's possible that
% there may still be problems
%
%
%
%% License:
%=====================================================================
%
% This is part of the Princeton MVPA toolbox, released under
% the GPL. See http://www.csbmb.princeton.edu/mvpa for more
% information.
% 
% The Princeton MVPA toolbox is available free and
% unsupported to those who might find it useful. We do not
% take any responsibility whatsoever for any problems that
% you have related to the use of the MVPA toolbox.
%
% ======================================================================
%
% NOTE: This function was written to allow for SPM5 compatability,
% and assumes SPM5 is installed and unmodified.  Specifically, this
% function makes use of .IMG and .HDR input/output functions in
% spm_dir/@nifti/private, strangely enough...

%% Defaults and setup
% single is used to set whether you'd like to open your data as a single
% precision instead of a double precision.  This allows you to save a
% signifigant amount of memory if you don't actually need double precision.
defaults.single = false;

%capture the arguements and populate the default values.
args = propval(varargin,defaults);


%% Check for spm functions
spmPres = which('spm_vol.m');
if isempty(spmPres)
  error('SPM not found.');
end

%% Mask Setup
maskvol = get_mat(subj,'mask',maskname);
mDims   = size(maskvol);%#ok<NASGU> %get the dimensions of the mask
mask    = find(maskvol);%get the relevant indexes of the mask (all non zero's)
mSize   = length(mask);%get the size of the mask

% check mask isn't empty
if isempty(mask)
  error('Empty mask passed to load_analyze_pattern()');
end


%% Initialize the data structure
subj = init_object(subj,'pattern',new_patname);

%% Convert filenames to a cell array
%if the file name is an array of characters
if ischar(filenames)
    
  if ~isempty(strfind(filenames,'*'))
    [pat,jnk,jnk] = fileparts(filenames); %#ok<NASGU>
    tmp = dir(filenames);
    filenames = {tmp(:).name};
    if ~isempty(pat)
      for i=1:length(filenames)
	filenames{i} = [pat '/' filenames{i}];
      end
    end   
  else
    filenames = {filenames};
  end
  
elseif ~iscell(filenames)
  error('Filenames are not in form of char or cell.');
end

nFiles = length(filenames);

%% Initialize the data structure
tmp_data = zeros(mSize ,nFiles); %#ok<NASGU>

disp(sprintf('Starting to load pattern from %i .IMG files',nFiles));

%% Create a volume structure
vol = spm_vol(filenames);
tmp_data = []; %#ok<NASGU>

    %%%%%%%%%%%%%%%%%%%%%%
    %sylvains contribution
    %%%%%%%%%%%%%%%%%%%%%%

total_m = 0;

% compute total number of EPI images
for h = 1:nFiles
  [m n] = size(vol{h}); %#ok<NASGU>
  total_m = total_m + m;
end;

% allocate all at once to avoid reshaping iteratively
tmp_data = zeros(mSize, total_m);
  
total_m = 0;

    %% end contribution

for h = 1:nFiles
  if mod(h,400)==0 %% added by AES because printout was annoying
    fprintf('\t%i',h);
  end
  
  [m n] = size(vol{h}); %#ok<NASGU>
  
  for i = 1:m
     curvol = vol{h}(i);
     
     % Enforce mask size
     if ~all(curvol.dim == size(maskvol))
       error('Inapropriate mask for .IMG');
     end

     % Load the data from the IMG file
     [Vdata] = spm_read_vols(curvol);
     
     % AES added option to global normalize
     if other_args.globalnorm
         Vdata=Vdata-other_args.globaltimecourse(h);
     end
  
     if args.single
       Vdata = single(Vdata);
     end
     
     tmp_subvol(1:mSize,i) = Vdata(mask);
     
  end
  
  % Reshape the data to be Voxels X Time
  %tmp_data(1:mSize,end+1:end+m) = tmp_subvol;
  
    %%%%%%%%%%%%%%%%%%%%%%
    %sylvains contribution
    %%%%%%%%%%%%%%%%%%%%%%
    tmp_data(1:mSize,total_m+1:total_m+m) = tmp_subvol;
    total_m = total_m + m;
    %% end contribution

end % for h

disp(' ');

%% Store the data in the pattern structure
subj = set_mat(subj,'pattern',new_patname,tmp_data);

%% Set the masked_by field in the pattern
subj = set_objfield(subj,'pattern',new_patname,'masked_by',maskname);

%% Add the history to the pattern
hist_str = sprintf('Pattern ''%s'' created by load_analyze_pattern',new_patname);
subj = add_history(subj,'pattern',new_patname,hist_str,true);

%% Add information to the new pattern's header, for future reference
subj = set_objsubfield(subj,'pattern',new_patname,'header', ...
			 'vol',vol,'ignore_absence',true);

%% Load the subject         
% This object was conceived under a tree. Store that information in
% the SUBJ structure
created.function = 'load_analyze_pattern';
created.args = args;
subj = add_created(subj,'pattern',new_patname,created);


end

