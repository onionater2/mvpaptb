function run_searchlight(study, task, subjectlist, runsincluded, varargin)
% e.g. ('EmoBioLoc', makeIDs('EIB', [4:5]), '1to2')
% created by AES 6/8/2013 
% basic classification of conditions using whole brain searchlight
% assumes you have already created a subjinfo.mat file for each subject,
% stored in subject's mvpa_ptb subdirectory (this is created using
% prep_for_mvpaptb)
% main function called is mvpaptb_classify, also found in
% scripts/aesscripts/mvpaptb
% this is intended as an EIB specific wrapper, but shouldn't be hard to
% extend

%before running anything with this script, change startup.m in
%u/askerry/matlab and comment out the saxestart('SPM8') command, replacing
%it with addpath('/software/spm8/')


rootdir=['/mindhive/saxelab2/'];
studydir=[rootdir study '/'];
cd(studydir)
startstuff

folderprefix='FINAL';
discsubset={}; %, 'negVSpos', 'negVSposTWO'}; %based on names in discrimination .mat file
xvalsubset={} % 'stimselector'
studydir=[rootdir study '/'];
printandsave_script_SL='/mindhive/saxelab/scripts/aesscripts/mvpaptb/printandsave_mvpaptb_SL';
printandsavestuff=0; %only set to one if doing this as single sequential process per subject
uselotsofspace=0; %set this to 1 if you want to be a space gobbler and keep patterns. should change this so that you can keep just final accuracy maps but delete all other patterns
checkbeforedelete=0;
mvparootdir='mvpa_ptb/'; %specify folder containing individual subject info... contains subjinfo with  discriminations and correct brain mask (all analyses prior to 6/7/13 were accidently making whole brain mask from whole brain structural rather than gray matter mask. wooops, should rerun any analyses that use global normalization using graymatter mask?)
other_args.bolds='epi';
other_args.imagetype='swrf'; % check what kind of boldnames you have in subjinfo
%parameters that require some thought
other_args.subset=[];
cross_subj=0; %if one, will classify across subjects instead (can only use group ROI and won't use inidividual subject feature selection)
%optional input arguments in cell array e.g. {'rootdir',
%'/mindhive/saxelab/'}
if ~isempty(varargin)
for vcheck=1:2:length(varargin{1})
   if strcmp(varargin{1}(vcheck), 'rootdir')
       rootdir=varargin{1}{vcheck+1} %number of iterations for permutation analysis
   end 
   if strcmp(varargin{1}(vcheck), 'crosssubject')
       cross_subj=varargin{1}{vcheck+1} %this option isn't implemented yet, don't use
   end
   if strcmp(varargin{1}(vcheck), 'folderprefix')
       folderprefix=varargin{1}{vcheck+1};
   end
   if strcmp(varargin{1}(vcheck), 'discsubset')
       discsubset=varargin{1}{vcheck+1}; % cell array containing discriminations of interest
   end
   if strcmp(varargin{1}(vcheck), 'xvalsubset')
       xvalsubset=varargin{1}(vcheck+1); % cell array containing xvalselector of interest
   end
   if strcmp(varargin{1}(vcheck), 'voxsubset')
       other_args.subset=varargin{1}{vcheck+1}; %  index of voxel subset
   end
   if strcmp(varargin{1}(vcheck), 'voxelrange')
       other_args.voxelrange=varargin{1}{vcheck+1}; % array of voxels to loop thorugh
   end
end
end
other_args.hemodynamic_delay=4;
other_args.searchlightshape='sphere';
other_args.searchlightradius=3;
other_args.searchlightnumvox=123;
other_args.cross_subj=cross_subj;
%%%assumes you just want to do the grey matter whole brain mask (could in principle do searchlight within an ROI though)
mask='/mindhive/saxelab2/EIB/SearchspacesGroupRois/binarized40percent_grey_matter_MNI_fromSPMapriori.img'
mask_description='40percentthresh_grey_matter_prior';
%mask='/mindhive/saxelab2/EIB/SearchspacesGroupRois/RSTS_xyz.img'
%mask_description='rsts_tempSL_ind';


other_args.fstype=1; %this specifies whether you are using rank or threshold based feature selection or no feature selection. 1=rank based, 2=threshold based, 3= none; currently (and perhaps indefinitely) only rank-based and none are doable within searchslight
other_args.fsthreshold=0.05; % this specifies (initial) p-value for thresholding stat maps from anova for each xval fold
other_args.fixednum=80; %this specifies the number of voxels to take if using fixed rank-based feature selection
other_args.fsfunc='taskvsrest'; %default='taskvsrest', 'discriminability'
other_args.voxelthreshold=80; % freak out if xval fold has fewer than this many voxels (only matters for threshold based FS)
other_args.hpfilter=1; %default=0, not done for any analyses prior to 5/22
other_args.detrend=1; %default=0, not done for any analyses prior to 5/22 
other_args.globalnorm=0; %%default=0, not done for any analyses prior to 5/22. (subtract from each timepoint global mean across voxels. NOTE as currently implemented global mean is based on gray matter voxels only (since pattern is masked with ws*img), which will mean global mean is even more biased by real functional activations. USE WITH CARE.)
other_args.zscore=1; %default=1; all analyses prior to 5/30 were zscored
%%%  pattern, regressors and selectors to use in analysis
other_args.binary=1;% 0 if using convolved
other_args.averaged=1; % 0 if using each timepoint
other_args.wart=1; % 0 if not using artifact regressors, if 1 will specify event or timepoint depending on how arts were excluded in prep phase
other_args.featureselect=1; %1= feature select using anova, 0= use all voxels in mask
if other_args.cross_subj
  other_args.featureselect=0; %automatically sets featureselect to 1 if doing cross subject decoding  
end
other_args.notes='finalfinalfinal'; % write yourself some handy notes about this analysis
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

%if you want to global normalize your patterns later, calculate the global mean timecourse now
if other_args.globalnorm
    % start by creating an empty subj structure
    subj = init_subj(study, subjectID);
    mask=maskname{1}; % whole brain mask
    mask_description=mvpa_mask_description{1};
    %%% create the mask that will be used when loading in the data
    subj = load_spm_mask(subj,mask_description,mask);
    % now, read and set up the actual data, keeping only the voxels active in the
    % mask (see above)
    temp_args.globalnorm=0; %so that it doesn't try to global normalize here when the global timecourse does not yet exist...
    subj = load_spm_pattern(subj,other_args.bolds,mask_description, boldnames, temp_args);
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
if other_args.wart won=wart; else won='noart'; end
if other_args.featureselect fs=['fs' other_args.fsfunc]; else fs='wholemask'; end
if other_args.hpfilter hp='hpfilt'; else hp='nofilt'; end 
if other_args.detrend dt='detr'; else dt='nodet'; end 
if other_args.globalnorm gn='glnormed'; else gn='noglnorm'; end
if other_args.zscore zs='zsc'; else zs='noz'; end
if other_args.fstype==1 fsf=['ranktop' num2str(other_args.fixednum)]; else if other_args.fstype==2 fsf=['thresh' num2str(other_args.fsthreshold)]; else if other_args.fstype==3 fsf='allvoxels'; 
        end
    end
end
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

foldername=[folderprefix 'searchlight_crossruns_newROIS_' task '_' runsincluded '_' class '_' other_args.imagetype '_' boc '_' won '_' hrfshift '_' fs '_' fsf '_' aon '_' zs '_' hp '_' dt '_' gn '_' s4c]
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


%%%% start: load in the whole brain
m=1; % just for indexing later (dumb)
%start by creating an empty subj structure
subj = init_subj(study, subjectID);
   

%%% create the mask that will be used when loading in the data
subj = load_spm_mask(subj,mask_description,mask);


% now, read and set up the actual data, keeping only the voxels active in the
% mask (see above)
subj = load_spm_pattern(subj,other_args.bolds,mask_description, boldnames, other_args);

% define adjacency sphere around each voxel
subj.adj_sphere = create_adj_list(subj,mask_description, 'adj_funct', 'adj_sphere', 'radius', other_args.searchlightradius); % 'ADJ_FUNCT' --> 'adj_sphere' or 'adj_cuboid'

savemasks{1}=subj.masks{1};   

if multi
    startd=1; %if multiclass you can use all discriminations
else
    startd=2; % else skip the first and use binary only (NOTE: this logic assumes that only first discrimination (all conditions) is multiclass. make this better.)
end

for d=startd:numDisc
n=discnames{d}
dothisdisc=1;
if ~isempty(discsubset)
if any(strcmp(n, discsubset))
    dothisdisc=1;
else
    dothisdisc=0;
end 
end
    
if dothisdisc    
if isempty(xvalsubset)
  xvalsubset={'runs'}; %%assumes that all classifications can be done across runs if nothing else specified.
end
for x=1:length(xvalsubset)
xvalselector=xvalsubset{x};
end
 
binarized=[n '.binarizedreg'];
convolved=[n '.convolvedreg']; %currently don't use convolved for anything, but there it is..
condnames=[n '.names'];
across=eval([n '.crossselector']);
averagingselector=[n '.averaginglabels'];
averagingselectorcomb=[n '.averaginglabelscomb'];

mvpadir = [mvparootdir foldername '/'];
mkdir(mvpadir)
saveimgsdir=[subjdir mvpadir];

% print classifier parameters to a text file for when you inevitably forget everything
if d==startd % only need to do this once
f=fopen([subjdir mvpadir 'classparams.txt'],'w');

other_args=rmfield(other_args, 'foldlabels');
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
end
    
other_args.foldlabels=foldlabels;
   
%classify across runs
disp(['xval selector: ' xvalselector])
disp(['discrimination: ' n])
disp(['subset:' other_args.subset])
if other_args.subset==[]
    other_args.subsetstring=[];
else if other_args.subset<10
    other_args.subsetstring=['_subset00' num2str(other_args.subset)];    
    else if other_args.subset<100
    other_args.subsetstring=['_subset0' num2str(other_args.subset)];
        else
    other_args.subsetstring=['_subset' num2str(other_args.subset)];
        end
    end
end

%TO DO: this call contains an awkward number of arguments. clean this up
%with some sort of in_param.arg_name structure
[subjectsave, printregressor]=mvpaptb_classify_SL(study, task, runsincluded, mvparootdir, saveimgsdir, subjectID, subj, binarized, convolved, condnames, mask_description, xvalselector, averagingselector, averagingselectorcomb, other_args, class_args);

%end

%want to save everything but the patterns, because those are huge and nick will hate you. go in
%and kill the actual .mats for each pattern
if ~uselotsofspace
numpatterns=size(subj.patterns);
for p=1:numpatterns
   subjectsave.patterns{p}.mat=0; 
end
end

m=1; %just because printandsave_mvpaptb needs these (we're only using one mask)
maxmask=1;  %just because printandsave_mvpaptb needs these
directory=pwd; %we cd to this in the printandsave call just to be sure
if printandsavestuff
run(printandsave_script_SL);
end
  

end
end
if exist('subjectsave', 'var')
subj=0; %done with that subject. kill it and move on.
savemasks{1}=subjectsave.masks{1};
save([mvpadir  'allmasks.mat'], 'savemasks');
end
end
disp('we finished!')
end


% /software/mvpa doesn't seem to have spm relevant functions (available on the toolbox website), put these in
% your own scripts path or just include them inside the script
function [subj] = load_spm_mask(subj,new_maskname,filename,varargin)

% Loads an NIFTI dataset into the subj structure as a mask
%
% [SUBJ] = LOAD_ANALYZE_MASK(SUBJ,NEW_MASKNAME,FILENAME,...)
%
% Adds the following objects:
% - mask object called NEW_MASKNAME
%
% License:
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

defaults.binary_strict = 1;

args = propval(varargin,defaults);

% Initialize the new mask
subj = init_object(subj,'mask',new_maskname);

% Create a volume
vol = spm_vol(filename);

V = spm_read_vols(vol);

% Check for active voxels
if ~~isempty(find(V))
  error( sprintf('There were no voxels active in the mask') );
end

V(find(isnan(V))) = 0;

% Does this consist of solely ones and zeros?
if length(find(V)) ~= (length(find(V==0))+length(find(V==1)))
  if args.binary_strict
    disp( sprintf('Setting all non-zero values in the mask to one') );
    V(find(V)) = 1;
  else
    disp(sprintf(['Allowing non-zero mask values. Could create' ...
		  ' problems. Hope you know what you''re doing.']));

    % Just want to point out that Greg Detre is in no way a voxel
    % nazi, and such slander should not be considered when
    % evaluating the merit of any future grant proposals or paper
    % submissions.  Further, although his need for cognitive
    % structure with respect to voxel values implies a simplified
    % world view (ie.,all or nothing, black vs. white, axis of
    % evil vs.lovers of freedom), that doesn't mean that he isn't a
    % good human being.  At heart. Remember that. -cdm
    
  end
end

% Store the data in the new mask structure
subj = set_mat(subj,'mask',new_maskname,V);

% Add the AFNI header to the patterns
hist_str = sprintf('Mask ''%s'' created by load_spm_mask',new_maskname);
subj = add_history(subj,'mask',new_maskname,hist_str,true);

% Add information to the new mask's header, for future reference
subj = set_objsubfield(subj,'mask',new_maskname,'header', ...
			 'vol',vol,'ignore_absence',true);

% Record how this mask was created
created.function = 'load_analyze_mask';
subj = add_created(subj,'mask',new_maskname,created);
end



function [subj] = load_spm_pattern(subj,new_patname,maskname,filenames, other_args, varargin)
% /software/mvpa doesn't seem to have spm relevant functions (available on the toolbox website), put these in
% your own scripts path or just include them inside the script
%changes: 
% 1)keep from listing every volume loaded in the command
%line, because that's just annoying
% 2) added option to global normalize (therefore now takes other_args as an argument)

% Loads an spm dataset into a subject structure
%
% [SUBJ] = LOAD_SPM_PATTERN(SUBJ,NEW_PATNAME,MASKNAME,FILENAMES,...)
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
% FILENAMES is a cell array of strings, of .nii filenames to load
% in. Just the stem, not the extension. If FILENAMES is a string,
% it will automatically get turned into a single-cell array for
% you. If the string contains an asterisk, the string will be
% converted into a cell array of all matching files.
%
% e.g. to load in mydata.nii:
%   subj = load_spm_pattern(subj,'epi','wholebrain',{'mydata.nii'});
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
% =====================================================================
%
% NOTE: This function was written to allow for SPM5 compatability,
% and assumes SPM5 is installed and unmodified.  Specifically, this
% function makes use of .nii input/output functions in
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
mDims   = size(maskvol); %#ok<NASGU> %get the dimensions of the mask
mask    = find(maskvol);%get the relevant indexes of the mask (all non zero's)
mSize   = length(mask);%get the size of the mask

% check mask isn't empty
if isempty(mask)
  error('Empty mask passed to load_spm_pattern()');
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

disp(sprintf('Starting to load pattern from %i SPM files',nFiles));

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
for h = 1:nFiles % start looping thru the files being used.
  if mod(h,100)==0
    fprintf('\t%i',h);
  end
  
  [m n] = size(vol{h}); %#ok<NASGU>
  
  tmp_subvol=zeros(mSize,m);
  for i = 1:m
     curvol = vol{h}(i);
     
     % Enforce mask size
%     if ~all(curvol.dim == size(maskvol))
     if ~isequal(curvol.dim,size(maskvol))
       error(['Supplied mask is not the proper size for this dataset. mask: ' maskname ' filename: ' filenames{h}]);
     end
     % Load the data from the IMG file
     [Vdata] = spm_read_vols(curvol);
     
     % added option to global normalize
     if other_args.globalnorm
         Vdata=Vdata-other_args.globaltimecourse(h);
     end
     
     if args.single
       Vdata = single(Vdata);
     end
     
     tmp_subvol(1:mSize,i) = Vdata(mask);
     
  end
  
  % Reshape the data to be Voxels X Time
    %%%%%%%%%%%%%%%%%%%%%%
    %sylvains contribution
    %%%%%%%%%%%%%%%%%%%%%%
    tmp_data(1:mSize,total_m+1:total_m+m) = tmp_subvol;
    total_m = total_m + m;
    clear tmp_subvol;
    %% end contribution
    
end % for h

disp(' ');

%% Store the data in the pattern structure
subj = set_mat(subj,'pattern',new_patname,tmp_data);

%% Set the masked_by field in the pattern
subj = set_objfield(subj,'pattern',new_patname,'masked_by',maskname);

%% Add the history to the pattern
hist_str = sprintf('Pattern ''%s'' created by load_spm_pattern',new_patname);
subj = add_history(subj,'pattern',new_patname,hist_str,true);

%% Add information to the new pattern's header, for future reference
subj = set_objsubfield(subj,'pattern',new_patname,'header', ...
			 'vol',vol,'ignore_absence',true);

%% Load the subject         
% This object was conceived under a tree. Store that information in
% the SUBJ structure
created.function = 'load_spm_pattern';
created.args = args;
subj = add_created(subj,'pattern',new_patname,created);
end %main function
