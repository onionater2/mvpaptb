function [scratch] = train_libsvm(trainpats,traintargs,in_args,cv_args)
% edited by AES to call propval_AES and convert trainpats to double as
% suggested here: https://groups.google.com/forum/?fromgroups=#!topic/mvpa-toolbox/VLnHXYC4zaI
% also edited by AES to add path to libsvm and exclude path to
% biolearning/svmtrain
% edited to add in in_args parameters to training_options

%AES: I believe this is equivalent to option 1 in liblinear: L2-regularized L2-loss support vector classification

% USAGE:
% [SCRATCH] = TRAIN_LIBSVM(TRAINPATS,TRAINTARGS,IN_ARGS,CV_ARGS)
%
% This is a support vector machine training function using the libsvm
% library. It train the classifier and makes is ready for testing.
%
% You need to call TEST_LIBSVM afterwards to assess how well this
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
%                             pass to svmtrain of libsm library.  See libsvm
%                             docs for more details. (example, '-s 0 - t 0 -c 1'
%                             for C-SVC with LINEAR kernel and a cost of 1).
%
%  in_args.search_for_c = (optional, default is 0).  Boolean, perform cross-validation
%                         within the training to find a -c value that gives best
%                         generalization within training set.  In this case, you cannot
%                         provide a -c or -v option in training_options and must also
%                         supply a k_fold_xval (see below).
%                    *****N.B. THIS OPTION NEEDS SOME WORKS. xxx
%
%  in_args.k_fold_xval = (must be supplied if search_for_c is true, default 0).  The number
%                        of k-folds to use for within-training_set cross-validation for
%                        optimizing c parameter.  See libsvm library doc for more info.
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
%
%======================================================================

% 02.27.09 - REHBM - created specific LIBSVM train function to
%                    distinguish from SVMLIGHT functions
%                    Based on train/test_svm.m that came with MVPA toolbox.
% 03.05.09 - REHBM - added support for c-parameter search using k-fold cross-validation

% TODO: args.training_options should be split into each possible argument with error checking?
%       k-fold cross validation should be done by splitting each run of the training set out. This
%          could be a recursive method by redefining the training set and training and testing.
%          But it requires that we extract run info from training targs.

%%
addpath('/mindhive/saxelab/scripts/aesscripts/svm/libsvm-3.17/matlab/') %% this is where the libsvm libary is located
rmpath('/software/matlab_versions/2010b/toolbox/bioinfo/biolearning') %% remove this directory from path because this contains another function called svmlearn

%% SORT ARGUMENTS
defaults.train_funct_name = 'train_libsvm';
defaults.test_funct_name  = 'test_libsvm';
defaults.training_options = '';
defaults.search_for_c     = 0;
defaults.k_fold_xval      = 0;
defaults.ignore_1ofn      = 'false';

% Args contains the default args, unless the user has over-ridden them
args = propval_AES(in_args,defaults);
scratch.class_args = args;
args = sanity_check(trainpats,traintargs,args);


%% Training Labels
%   Labels will be 1:K for K conditions
[train_max_val trainlabs] = max(traintargs); % trainlabs = training labels, max index of max val


%% Parameter search?
% should we search for the best c-parameter by cross-
% validation within the training set.
% N.B. This probbaly needs to account for the non-
%      independence across consecutive TRs, which are
%      also most likely of the same class for a block
%      design.

if args.search_for_c
    % make sure user didn't try to set a specific c with -c in the args.training_options
    if ~isempty(regexp(args.training_options,'.*-[cv].*','once'))
       error('train_libsvm.m - cannot supply -c or -v option when search_for_c is true (1)')
    end
    if args.k_fold_xval < 1
        error('train_libsvm.m - k_fold_xval (%d) must be supplied when search_for_c is true (1)',args.k_fold_xval)
    end
    starttime = tic;
    csearch.c         = 2.^(-7:1);
    csearch.accuracy  = repmat(NaN,[1 length(csearch.c)]);
    csearch.loop_time = repmat(NaN,[1 length(csearch.c)]);    
    for i = 1:length(csearch.c)       
        looptime = tic;
        c = num2str(csearch.c(i));
        argstring = [args.training_options ' -v ' num2str(args.k_fold_xval) ' -c ' c];
        csearch.accuracy(i) = svmtrain(trainlabs',trainpats',argstring);
        csearch.loop_time(i) = toc(looptime);
    end
    csearch.total_time = toc(starttime);

    % get c with best performance
    csearch.best_c = csearch.c(find(csearch.accuracy(1,:)==max(csearch.accuracy(1,:)),1));

    % update args.training_options
    %args.training_options = [args.training_options ' -c ' num2str(csearch.best_c)];
    args.training_options = [args.training_options ' -c ' num2str(csearch.best_c) ' -s ' num2str(in_args.svm_type) ' -t ' num2str(in_args.kernel_type)]; % AES edit
else
    args.training_options = [args.training_options ' -c ' num2str(in_args.cost) ' -s ' num2str(in_args.svm_type) ' -t ' num2str(in_args.kernel_type)]; %AES edit
end

args.training_options

%% Train the classifier
[scratch.model] = svmtrain(trainlabs',sparse(double(trainpats))',args.training_options);


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
