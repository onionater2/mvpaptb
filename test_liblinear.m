function [acts scratchpad] = test_liblinear(testpats,testtargs,scratchpad) 

% Testing function for LIBLINEAR (linear svm) 
% 
% [ACTS SCRATCH] = TEST_LIBLINEAR(TESTPATS,TESTTARGS,SCRATCH) 
% 
% This is the testing function that fits with TRAIN_LIBLINEAR. See that 
% for more info 
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
%====================================================================== 

% 02.27.09 - REHBM - created specific LIBSVM train function to 
%                    distinguish from SVMLIGHT functions 
%                    Based on train/test_svm.m that came with MVPA toolbox. 


%% validate arguments 
if ~exist('scratchpad','var') 
  scratchpad = []; 
end 

sanity_check(testpats,testtargs,scratchpad); 


%% Test Labels 
%   Labels will be 1:K for K conditions 
[test_max_val testlabs]  = max(testtargs); % testlabs = testing labels, max index of max val 

% Now test generalization performance on the test data 
% Classify the test set using svmpredict 
[scratchpad.predicted_label scratchpad.accuracy scratchpad.decision_values] = ... 
    predict(testlabs', sparse(testpats)', scratchpad.model, '-b 0'); 


%% trying to make LIBLINEAR compatible with default perfmat functions. 
% 
% ACTS = is an nOuts x nTestTimepoints matrix that contains the activations of the output units at test. 
% We'll fill acts using winner-take-all activation based on scratchpad.predicted_labels 
acts = zeros(size(testtargs)); % initilize to all zeros 
for i = 1:size(acts,1) 
    acts(i,scratchpad.predicted_label==i) = 1; % otherwise it remains zero from initilization 
end 


%% sanity_check 
function [] = sanity_check(testpats,testtargs,scratchpad) 


if size(testpats,2) ~= size(testtargs,2) 
  error('Different number of testing pats and targs timepoints'); 
end 
