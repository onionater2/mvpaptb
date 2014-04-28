function [scratch] = train_liblinear (trainpats,traintargs,in_args,cv_args) 
% edited by AES to call propval_AES and convert trainpats to double as
% suggested here: https://groups.google.com/forum/?fromgroups=#!topic/mvpa-toolbox/VLnHXYC4zaI
% also edited by AES to add path to liblinear and exclude path to
% biolearning/svmtrain
% edited to add in in_args parameters to training_options

% use the following key to assign values in run_classification_EIB.m
% 	 0 -- L2-regularized logistic regression (primal)\n"
% 	 1 -- L2-regularized L2-loss support vector classification (dual)\n"	
% 	 2 -- L2-regularized L2-loss support vector classification (primal)\n"
% 	 3 -- L2-regularized L1-loss support vector classification (dual)\n"
% 	 4 -- support vector classification by Crammer and Singer\n"
% 	 5 -- L1-regularized L2-loss support vector classification\n"
% 	 6 -- L1-regularized logistic regression\n"
% 	 7 -- L2-regularized logistic regression (dual)or regression\n"
% 	11 -- L2-regularized L2-loss support vector regression (primal)\n"
% 	12 -- L2-regularized L2-loss support vector regression (dual)\n"
% 	13 -- L2-regularized L1-loss support vector regression (dual)\n"
    % note to self: as per Andrew Ng's Feature selection, l1 vs l2 regularization, and rotational invariance paper, expect l1 regularization be better than l2 regularization if you have a lot less examples than features.
    

% USAGE : 
% [SCRATCH] = TRAIN_LIBLINEAR(TRAINPATS,TRAINTARGS,IN_ARGS,CV_ARGS) 
% 
% This is a linear support vector machine training function using the 
% liblinear library. It train the classifier and makes is ready for testing. 
% 
% You need to call TEST_LIBLINEAR afterwards to assess how well this 
% generalizes to the test data. 


% PATS = nFeatures x nTimepoints 
% TARGS = nOuts x nTimepoints 
% 
% SCRATCH contains all the other information that you might need when 
% analysing the network's output, most of which is specific to 
% backprop. Some of this information is redundantly stored in multiple 
% places. This gets referred to as SCRATCHPAD outside this function 
% 
% The classifier functions use a IN_ARGS structure to store possible 
% arguments (rather than a varargin and property/value pairs). This 
% tends to be easier to manage when lots of arguments are 
% involved. 
% 
% IN_ARGS are the various arguments that can be passed in for type 
% of kernels used and the learning parameters. 
% 
% IN_ARGS: 
%  in_args.training_options = (optional, default is '').  String of options to 
%                             pass to svmtrain of liblinear library. See liblinear 
%                             docs for more details. (example, '-s 1 - c 1' for 
%                             L2-loss SVM with LINEAR kernel and a cost of 1). 
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

% 03.01.09 - REHBM - created specific LIBLINEAR train function to 
%                    distinguish from SVMLIGHT/LIBSVM functions 
%                    Based on train/test_svm.m that came with MVPA toolbox. 

addpath('/mindhive/saxelab/scripts/aesscripts/svm/liblinear-1.93/matlab/') %% this is where the liblinear libary is located
rmpath('/software/matlab_versions/2010b/toolbox/bioinfo/biolearning') %% remove this directory from path because this contains another function called svmlearn



%% validate argument 
defaults.train_funct_name = 'train_liblinear'; 
defaults.test_funct_name  = 'test_liblinear'; 
defaults.training_options = ''; 
defaults.ignore_1ofn      = 'false'; 

% Args contains the default args, unless the user has over-ridden them 
args = propval_AES(in_args,defaults); 
scratch.class_args = args; 
args = sanity_check(trainpats,traintargs,args); 



%% Training Labels 
%   Labels will be 1:K for K conditions 
[train_max_val trainlabs] = max(traintargs); % trainlabs = training labels, max index of max val 

args.training_options = [args.training_options ' -c ' num2str(in_args.cost) ' -s ' num2str(in_args.svm_type)]; %AES edit
args.training_options

%% *** TRAINING THE CLASSIFIER... *** 
[scratch.model] = train(trainlabs',sparse(double(trainpats))',args.training_options);


%% Local Functions 
function [args] = sanity_check(trainpats,traintargs,args) 

if size(trainpats,2)==1 
  error('Can''t classify a single timepoint'); 
end 

if size(trainpats,2) ~= size(traintargs,2) 
  error('Different number of training pats and targs timepoints'); 
end 

[isbool isrest isoveractive] = check_1ofn_regressors(traintargs); 
if ~isbool || isrest || isoveractive 
  if ~args.ignore_1ofn 
    warning('Not 1-of-n regressors'); 
  end 
end 