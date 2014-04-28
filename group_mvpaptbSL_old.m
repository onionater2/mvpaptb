function group_mvpaptbSL(subjectlist, classfolder, task, runsincluded, varargin)
%this script is a disaster and needs lots of little in script changes for
%the particular contrast you are running. FIX IT! (in the mean time, search "change this" to find all lines that need specific updating)
%created by AES on 6/27 for group analysis of searchlight accuracy mapes
%within subjects averages across folds and creates crossfold avg map in
%each subject's directory
%then averages across subjects and makes group stat map (along with avg acc
%map and std map
%optional argument: {'chance', 0.25}
numSubj=length(subjectlist);
firstsubj=subjectlist{1}
lastsubj=subjectlist{numSubj}
subjectrange=['subj_' firstsubj(end-1:end) 'to' lastsubj(end-1:end)];
runRFX=1;

studydir='/mindhive/saxelab2/EIB/';
groupdir=([studydir '/RandomEffects/group_' classfolder, subjectrange '/'])
maskfile='/mindhive/saxelab2/EIB/SearchspacesGroupRois/binarized40percent_grey_matter_MNI_fromSPMapriori.img'


contrastdetails(1).name='negfVSposf';
contrastdetails(2).name='negcVSposc';
contrastdetails(3).name='negVSposONE';
contrastdetails(4).name='negVSposTWO';
contrastdetails(1).selectors={'runs', 'crossmatchedruns'};
contrastdetails(2).selectors={'runs', 'crossmatchedruns'};
contrastdetails(3).selectors={'other'};
contrastdetails(4).selectors={'other'};
contrastdetails(1).selectors={'runs', 'crossmatchedruns'};
contrastdetails(2).selectors={'runs', 'crossmatchedruns'};
contrastdetails(3).selectors={'other'};
contrastdetails(4).selectors={'other'};


mkdir(groupdir)
nonparametric=0; % see beginning of nonparametric impliementation below (but signrank is weird in matlab...)

chance=0.5; %default
if size(varargin)>0
    if strcmp(varargin{1}, 'chance')
        chance=varargin{2}
    else
        error('unknown parameter as argument')
    end
end

%just for subject 1, figure out the relevant discriminations (this assumes
%all subjects have same discriminations
       subjectID=subjectlist{1};
       subjdir=[studydir subjectID '/'];
       cd(subjdir)
       disc=load([subjdir 'mvpa_ptb/discriminations_' task '_' runsincluded]);
       discnames=fieldnames(disc)
       numDisc=length(discnames);
       selector{1}='runs';
       % selector{2}='crossmatchedruns'; %change this
       selector{2}='selector';
%go through each discrimination
       for d=8%10%8%7%1:numDisc %%change this d=negfposf, 
       discname=discnames{d};
            for c=1:2
            if c==1 
                selectorname='runs'; 
                realselectornames{1}='runs';
                %realselectornames{1}='crossmatchedruns'; %change this
            else
                selectorname='other';
                selectorcount=1;
                tempimgs=dir([subjdir 'mvpa_ptb/' classfolder '/' discname '*' selector{c} '*.img']);
                for t=1:length(tempimgs)
                    name=tempimgs(t).name;
                    findtrain=findstr(name, 'train');
                    trainvar=name(findtrain+5: findtrain+8);
                    found=0;
                    for x=1:length(realselectornames)
                    if strcmp(realselectornames{x},trainvar)
                    found=trainvar;
                    end
                    end
                    if ~found
                    selectorcount=selectorcount+1;
                    realselectornames{selectorcount}=trainvar;
                    end
                end
            end
            end
            for c=2:length(realselectornames) %%% change this
                realselectorname=realselectornames{c};
           % filestring=[discname '*' realselectorname '*.img'];
            %filestring='negfVSposf.binarizedreg_runs_trainruns_srchacc_*.img'
            %filestring='negcVSposc.binarizedreg_runs_trainruns_srchacc_*.img'
            %filestring='negfVSposf.binarizedreg_crossmatchedruns_trainruns_srchacc_*.img'
            filestring=['negVSposTWO.binarizedreg_crossrunsTWOselector_tra*' realselectorname '*.img'] %%change this

%load mask
maskfileinfo=spm_vol(maskfile);
mask=spm_read_vols(maskfileinfo);            

%go through each subject and make their crossfold avg .imgs
subjcount=0;
for s=1:length(subjectlist)
       subjectID=subjectlist{s};
       subjdir=[studydir subjectID '/'];
       cd(subjdir)
            mvpadir=[subjdir 'mvpa_ptb/' classfolder '/'];
            cd(mvpadir)
            imgnames=dir(filestring);
            numfolds=size(imgnames,1)/2;
            if numfolds>0
                subjcount=subjcount+1;
            imgstring=imgnames(1).name(1:end-11);
            p=spm_select('list', mvpadir, ['^' imgstring '.*\.img$']); %% not sure why you need this special filter instead of swrf*.img, but you do

                files = spm_vol(p);
                for i=1:numfolds
                   i
                   files(i).fname
                if files(i).fname(end-19:end-14) ~='minus0'
                data(:,:,:,i)=spm_read_vols(files(i));
                data(data==0)=NaN;
                end
                end
                voxelwiseAvg(:,:,:,subjcount)=mean(data,4).*mask;; % 
                voxelwiseDiffFromChance(:,:,:,subjcount)=(mean(data,4)-chance).*mask;
                clearvars data
                
                
            writeTemplate=files(1); % getting template from a single fold classification .img, keeping pinfo
            classificationpinfo=writeTemplate.pinfo
            writeTemplate.dt = [spm_type('float64') spm_platform('bigend')];
            writeTemplate.fname = [subjectID '_' discname '_' realselectorname '_crossfoldMEANACC.img'];
            avg_singleoutput.fname = [subjectID '_' discname '_' realselectorname '_crossfoldMEANACC.img'];
            avg_singleoutput=spm_create_vol(writeTemplate);
            avg_singleoutput = spm_write_vol(avg_singleoutput, voxelwiseAvg(:,:,:,subjcount));        
            copyfile([subjectID '_' discname '_' realselectorname '_crossfoldMEANACC.img'], [groupdir 'IND_' subjectID '_' discname '_' realselectorname '_crossfoldMEANACC.img'])
            copyfile([subjectID '_' discname '_' realselectorname '_crossfoldMEANACC.hdr'], [groupdir 'IND_' subjectID '_' discname '_' realselectorname '_crossfoldMEANACC.hdr'])
           
            writeTemplate.fname = [subjectID '_' discname '_' realselectorname '_crossfoldDiffChance.img'];
            avg_singleoutput.fname = [subjectID '_' discname '_' realselectorname '_crossfoldDiffChance.img'];
            avg_singleoutput=spm_create_vol(writeTemplate);
            avg_singleoutput = spm_write_vol(avg_singleoutput, voxelwiseDiffFromChance(:,:,:,subjcount));        
            copyfile([subjectID '_' discname '_' realselectorname '_crossfoldDiffChance.img'], [groupdir 'IND_' subjectID '_' discname '_' realselectorname '_crossfoldDiffChance.img'])
            copyfile([subjectID '_' discname '_' realselectorname '_crossfoldDiffChance.hdr'], [groupdir 'IND_' subjectID '_' discname '_' realselectorname '_crossfoldDiffChance.hdr'])
 
            end
            
end

%now do group analysis if there is at least one subject
if runRFX
if subjcount>1;
cd(groupdir)

if nonparametric
dim=size(voxelwiseAvg)
groupW=zeros(dim(1:3));
 for x=1:dim(1)
     for y=1:dim(2)
         for z=1:dim(3)
             for subj=1:dim(4)
                 vector(subj)=voxelwiseAvg(x,y,z,subj);
             end
             if sum(isnan(vector))==0
             [p, h, stat]=signrank(vector, .5);
             stat.signedrank
             groupW(x,y,z)=stat.signedrank;
             end
         end
     end
 end
end
                groupAvg=nanmean(voxelwiseAvg,4);
                groupSTD=nanstd(voxelwiseDiffFromChance,0,4);
                groupSTD(groupSTD==0)=NaN;
                groupAvg(groupAvg==0)=NaN;
                groupSE=groupSTD./sqrt(subjcount);
                maskgroup=mask(:,:,:,1);

                
                groupSE=groupSE.*mask;
                groupDiff=groupAvg-chance;
                groupT=groupDiff./groupSE;
                groupDiff=groupDiff.*mask;
                groupT=groupT.*mask;
                
            inputfiles={'groupAvg', 'groupDiff', 'groupT', 'groupSE',}; 
            outputfiles={'AVG_output', 'Diff_output', 'T_output', 'SE_output',};  
            outputfilenames={'_crosssubjMEANACC.img','_crosssubjACCminus50.img', '_crosssubjT.img', '_crosssubjSE.img',};
            for x=1:length(outputfiles)
            outputfile=outputfiles{x};    
            writeTemplate.fname = [discname '_' realselectorname outputfilenames{x}];
            name=[discname '_', realselectorname outputfilenames{x}];
            eval([outputfile '.fname = name;'])
            eval([outputfile, '=spm_create_vol(writeTemplate);'])
            eval([outputfile,' = spm_write_vol(' outputfile ', ',inputfiles{x},'(:,:,:));']) 
            
            end
end
end
            end
            end
end

