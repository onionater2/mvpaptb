function SL_AvgAcrossPartitions(subjectlist, classfolder, task, runsincluded)
studydir='/mindhive/saxelab2/EIB/'
numSubj=length(subjectlist);
firstsubj=subjectlist{1}
lastsubj=subjectlist{numSubj}
subjectrange=['subj_' firstsubj(end-1:end) 'to' lastsubj(end-1:end)];
groupdir=([studydir '/RandomEffects/group_' classfolder, subjectrange '/'])
maskfile='/mindhive/saxelab2/EIB/SearchspacesGroupRois/binarized40percent_grey_matter_MNI_fromSPMapriori.img';
%load mask
maskfileinfo=spm_vol(maskfile);
mask=spm_read_vols(maskfileinfo);
%%
cd(groupdir)
contypes={'DiffChance', 'MEANACC'}
for s=1:numSubj
    subjectID=subjectlist{s};
for c=1:length(contypes)
    contype=contypes{c}
    subjfiles=dir(['IND_' subjectID '*crossruns*selector*crossfold' contype '.img'])
    for f=1:length(subjfiles)
       p=subjfiles(f).name
       file = spm_vol(p); 
       data(:,:,:,f)=spm_read_vols(file);
       data(data==0)=NaN;     
    end
    voxelwiseAvg(:,:,:)=mean(data,4).*mask;

            writeTemplate=file; % getting template from a single partition classification .img, keeping pinfo
            classificationpinfo=writeTemplate.pinfo
            writeTemplate.dt = [spm_type('float64') spm_platform('bigend')];
            writeTemplate.fname = ['IND_' subjectID '_negVSposABSTRACT_partitionsAvgd_crossfold' contype '.img'];
            avg_singleoutput.fname = ['IND_' subjectID '_negVSposABSTRACT_partitionsAvgd_crossfold' contype '.img'];
            avg_singleoutput=spm_create_vol(writeTemplate);
            avg_singleoutput = spm_write_vol(avg_singleoutput, voxelwiseAvg(:,:,:));  
end
end

