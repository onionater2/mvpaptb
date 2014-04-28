function addclassparamslater(subjectlist)

%analysisfolder='gnbgnb0_swrf_binary_wart_featureselect_averaged';
analysisfolder='libsvmlibsvm0_swrf_binary_wart_featureselect_averaged_costoptimized'
other_args.bolds='epi';
other_args.imagetype='swrf'; % check what kind of boldnames you have in subjinfo
other_args.fsthreshold=0.05; % this specifies (initial) p-value for thresholding stat maps from anova for each xval fold
other_args.fsfunc='anova';
other_args.voxelthreshold=20; % freak out if xval fold has fewer than this many voxels
other_args.binary=1;% 0 if using convolved
other_args.averaged=1; % 0 if using each timepoint
other_args.wart=1; % 0 if not using artifact regressors
other_args.featureselect=1; %1= feature select using anova, 0= use all voxels in mask
other_args.notes='ridge regression (I think L2-regularized L1-loss)'; % special notes to self about this analysis
other_args.classifier='libsvm' % 'lda' 'ridge', 'gnb'
    %class_args.nHidden = 0; %specifies number of hidden layers
    %class_args.search4c=2; % just specifies that no cost parameter is used here (optimized or fixed)
    %class_args.penalty = 0.05; %specificies how to penalize low weights. will multiply this constant times the number of voxels in the last xval fold the more voxels the more you want to penalize   
    class_args.search4c=1;
    class_args.kernel_type = 0; %0 -- linear: u'*v; 1 -- polynomial: (gamma*u'*v + coef0)^degree; 2 -- radial basis function: exp(-gamma*|u-v|^2); 3 -- sigmoid: tanh(gamma*u'*v + coef0)  
    %class_args.cost = 1;
    %class_args.search_for_c = 1; % 1= search for optimal c by xvalidating within training set, 0 = use .cost
    %class_args.search4c=1; %  specifies that  cost parameter is optimized within the training set
    class_args.k_fold_xval      = 8; %how many xval folds to do within training set to find optimal C parameter
    class_args.svm_type = '0 (C-SVC)'; %svm_type : set type of SVM (default 0): 0 -- C-SVC; 1 -- nu-SVC; 2 -- one-class SVM; 3 -- epsilon-SVR; 4 -- nu-SVR
    %class_args.cost = 1; %unclear why liblinear doesn't have search_for_c option like libsvm does



% analysisfolder='bpbp0_swrf_binary_wart_featureselect_averaged';
% 
% other_args.bolds='epi';
% other_args.imagetype='swrf'; % check what kind of boldnames you have in subjinfo
% other_args.fsthreshold=0.05; % this specifies (initial) p-value for thresholding stat maps from anova for each xval fold
% other_args.fsfunc='anova';
% other_args.voxelthreshold=20; % freak out if xval fold has fewer than this many voxels
% other_args.binary=1;% 0 if using convolved
% other_args.averaged=1; % 0 if using each timepoint
% other_args.wart=1; % 0 if not using artifact regressors
% other_args.featureselect=1; %1= feature select using anova, 0= use all voxels in mask
% other_args.notes='need to figure out if 2 layer with backprop is linear. what is the activation function?'; % special notes to self about this analysis
% other_args.classifier='bp' % 'lda' 'ridge', 'gnb'
%     class_args.nHidden = 0; %specifies number of hidden layers
%     %class_args.search4c=2; % just specifies that no cost parameter is used here (optimized or fixed)
%     %class_args.penalty = 0.05; %specificies how to penalize low weights. will multiply this constant times the number of voxels in the last xval fold the more voxels the more you want to penalize   
%     %class_args.search4c=0;
%     %class_args.kernel_type = 0; %0 -- linear: u'*v; 1 -- polynomial: (gamma*u'*v + coef0)^degree; 2 -- radial basis function: exp(-gamma*|u-v|^2); 3 -- sigmoid: tanh(gamma*u'*v + coef0)  
%     %class_args.cost = 1;
%     %class_args.search_for_c = 1; % 1= search for optimal c by xvalidating within training set, 0 = use .cost
%     %class_args.search4c=1; %  specifies that  cost parameter is optimized within the training set
%     %class_args.k_fold_xval      = 8; %how many xval folds to do within training set to find optimal C parameter
%     %class_args.svm_type = 0; %svm_type : set type of SVM (default 0): 0 -- C-SVC; 1 -- nu-SVC; 2 -- one-class SVM; 3 -- epsilon-SVR; 4 -- nu-SVR
%     %class_args.cost = 1; %unclear why liblinear doesn't have search_for_c option like libsvm does
% 




for s=1:length(subjectlist)
    
    subject=subjectlist{s};
    
subjdir=['/mindhive/saxelab2/EIB/' subject '/mvpa_ptb/' analysisfolder '/']

% print classifier parameters to a text file
f=fopen([subjdir 'classparams.txt'],'w');

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
end
