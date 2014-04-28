function run_grid_searchlight(subjectlist,doitinthismanychunks,numvoxelsInSL)
%e.g. run_grid_searchlight(makeIDs('EIB', [1]), 20, 163962)
%RSTS_xyz has 3377 voxels (for debugging)
%binarzed 40% has 163962





discriminations={'negfVSposf'}; %{'negVSposONE', 'negVSposTWO', 'negfVSposf', 'negcVSposc'}
xvalsubsets={'runs'};%'crossmatchedruns'};%{'crossrunsONEselector', 'crossrunsTWOselector', 'runs', 'runs'}

chunksize=floor(numvoxelsInSL/doitinthismanychunks);
startvox=1;
for i=1:doitinthismanychunks-1
thisrange=[startvox:startvox+chunksize];
startvox=startvox+chunksize+1;
voxelranges{i}=thisrange;
end
voxelranges{doitinthismanychunks}=[startvox:numvoxelsInSL];

% cd(['/mindhive/saxelab2/EIB/' 'SAX_EIB_13' '/mvpa_ptb/FINALsearchlight_crossruns_newROIS_EIB_main_1to8_libsvm_swrf_bin_tpoint_hrfshift4_fsactivity_ranktop80_avgd_zsc_hpfilt_detr_noglnorm_costfix']) 
% waititout=1
% while waititout
%    checkit=dir('negcVSposc.binarizedreg_crossmatchedruns_trainruns_srchacc_1_subset024.img');
% %checkit=dir('negVSposONE.binarizedreg_crossrunsONEselector_trainoddpersons_srchacc_1_subset029.img');
%    if length(checkit)>0
%        waititout=0;
%    end
%    pause(60)
% end

for subj=1:length(subjectlist)
    subject=subjectlist{subj};
    subjnum=str2num(subject(end-1:end));
    subjnum=subjnum-1;
    lastsubject=['SAX_EIB_' num2str(subjnum)]
 %   if subj>1
% cd(['/mindhive/saxelab2/EIB/' 'SAX_EIB_19' '/mvpa_ptb/FINALsearchlight_crossruns_newROIS_EIB_main_1to8_libsvm_swrf_bin_tpoint_hrfshift4_fsactivity_ranktop80_avgd_zsc_hpfilt_detr_noglnorm_costfix']) 
% waititout=1
% while waititout
%    checkit=dir('negfVSposf.binarizedreg_runs_trainruns_srchacc_1_subset024.img');
% %checkit=dir('negVSposONE.binarizedreg_crossrunsONEselector_trainoddpersons_srchacc_1_subset029.img');
%    if length(checkit)>0
%        waititout=0;
%    end
%    pause(60)
% end
 %   end
    
for d=1:length(discriminations)
    disc=discriminations{d}
    xval=xvalsubsets{d}    
    for iter=1:doitinthismanychunks 
    voxelrange=voxelranges{iter}
command=['run_searchlight(''EIB'', ''EIB_main'', {''' subject '''}, ''1to8'', {''discsubset'', ''' disc ''', ''xvalsubset'', ''' xval ''', ''voxsubset'', ' num2str(iter) ', ''voxelrange'', [' num2str(voxelrange) ']})']   
name= ['SL_' num2str(iter) '_' subject]
gridSubmitAES(command, name)
    end
end
end
end