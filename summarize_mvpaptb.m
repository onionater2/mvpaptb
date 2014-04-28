function summarize_mvpaptb(study, subjectlist, mvparesults, task, runsincluded)
%created by AES on 5/1/2013 to generate hasty (read ugly) summary of an
%mvpaptb classification output.

numSubj=length(subjectlist);
pdim=ceil(sqrt(numSubj));
%indMVPAdir='mvpa_ptb_old/';
indMVPAdir='mvpa_ptb/';
localizertasks={'tomloc', 'EmoBioLoc', 'group'}

firstsubj=subjectlist{1}
lastsubj=subjectlist{numSubj}
subjectrange=['_subj' firstsubj(end-1:end) 'to' lastsubj(end-1:end)];
minsubj2analyze=6;

rootdir=['/mindhive/saxelab2/' study '/'];
cd(rootdir)
mkdir([rootdir 'mvpaptb/' mvparesults, subjectrange 'xx/']);
mkdir([rootdir 'mvpaptb/' mvparesults, subjectrange, 'xx/specifics/']);
subjmvpa=[indMVPAdir mvparesults '/'];
groupmvpadir=[rootdir 'mvpaptb/' mvparesults, subjectrange 'xx/'];

%roilist={'lSTS_peelenpeak_xyz'};
%roilist={'rFFA_kanparcel_EmoBioLoc',...
%    'rpSTS_BDbiomot',...
%    'rSTS_kanparcel_EmoBioLoc',...
%    'rOFA_kanparcel_EmoBioLoc'}
%roilist={'whole_brain_skullstripped'};
roilist={
    'RTPJ_tomloc',...
    'LTPJ_tomloc',...
    'RSTS_tomloc',...
    'LSTS_tomloc',...
    'PC_tomloc',...
    'rFFA_kanparcel_EmoBioLoc',...
    'rpSTS_BDbiomot',...
    'rSTS_kanparcel_EmoBioLoc',...
    'rOFA_kanparcel_EmoBioLoc',...
    'MMPFC_tomloc',...
    'DMPFC_tomloc',...
    'VMPFC_tomloc',...
    'lSTS_peelenpeak_xyz',...
    'rSTS_peelenflip_xyz'};
    %'lSTS_kanparcel_EmoBioLoc',...
    %'lFFA_kanparcel_EmoBioLoc',...
    %'MPFC_peelenpeak_xyz',...
    %'lOFA_kanparcel_EmoBioLoc'
%     'rLOC_found_EmoBioLoc',...
%     'lLOC_found_EmoBioLoc',...
%     'rinsula_wfu_xyz',...
%     'rvSTR_reward_xyz',...
%     'vmPFC_reward_xyz',...
%     'right_ant_temporal_xyz',...
%     'ramygdala_wfu_xyz',...
%     'lvSTR_reward_xyz',... 
%     'linsula_wfu_xyz',...
%     'left_ant_temporal_xyz',...
%     'lamygdala_wfu_xyz'


% roilist={
%     'DM11_group.img',...
%     'DM21_group.img',...
%     'DM12_group.img',...
%     'DM22_group.img',...
%     'DV11_group.img',...
%     'DV21_group.img',...
%     'DV12_group.img',...
%     'DV22_group.img',...
%     'MV11_group.img',...
%     'MV21_group.img',...
%     'MV12_group.img',...
%     'MV22_group.img',...
%     'DMPFC_group.img',...
%     'MMPFC_group.img',...
%     'VMPFC_group.img',...
%     };



numROIs=length(roilist);

% figure out possible discriminations
template=load([rootdir subjectlist{1} '/' indMVPAdir 'discriminations_' task '_' runsincluded]);
templatesubjinfo=load([rootdir subjectlist{1} '/' indMVPAdir 'subjinfo_' task '_' runsincluded]);
condarray=templatesubjinfo.allcond_names
condstring=condarray{1};
for x=2:length(condarray) 
    condstring=[condstring '_' condarray{x}];
end

ds=fieldnames(template)
numDisc=length(ds);
        for d=1:numDisc
            desc=ds{d};

            % figure out if this disc has multiple selector options
            thedir=[rootdir subjectlist{1} '/' subjmvpa roilist{1} '*' desc '*classification.mat'];
            templateDES=dir(thedir); % to identify num of disc variants for the descrimination
            numDES=size(templateDES);
            for dse=1:numDES(1)
            tempdesc=templateDES(dse).name;
            roicount=0;
            hasroi=zeros(numROIs, length(subjectlist));
            for roin=1:numROIs
                roiname=roilist{roin};
               subjcount=0; 
               index=strfind(tempdesc, desc);
               tempdesc=tempdesc(index:end);
               disp(['working on discrimination ' tempdesc ' in roi ' roilist{roin}]);
               for subj=1:numSubj
                subject=subjectlist{subj};
                subjnum=(subject(end-1:end));
                subjdir=[rootdir subject '/'];
                resultsdir=[subjdir subjmvpa];
                cd(resultsdir)
                findfolder=dir(resultsdir);
                if isempty(findfolder)
                    disp(['couldn"t find folder for subject ' subject])
                else
                  
                    subjdesc=dir([resultsdir roiname '*' tempdesc]);
                    if isempty(subjdesc)
                        disp(['subject ' subject ' doesnt have that discrimination for roi ' roiname])
                        hasroi(roin,subj)=0;
                    else
                    currentdesc=subjdesc.name;
                    currentdesc=currentdesc(1:end-19);
                        if sum(hasroi(roin,:))==0
                           roicount=roicount+1;
                           roinames{roicount}=roiname;
                        end
                        hasroi(roin,subj)=1;
                        subjcount=subjcount+1;
                        subjlabels{roin}{subjcount}=subjnum; 

              %% do all the relevant things  
                
              %% get numvoxels
%%%%              load([currentdesc 'voxelsINmasks.mat']);
%%%%              voxelsGroup{subjcount}=masksize;
              %% get vector of classification accuracies across folds
              load([currentdesc '_classification.mat']);
              numIterations=length(results.iterations);
              numConds=length(unique(results.iterations(1).perfmet.desireds));
              chance=1/numConds;
              for it=1:numIterations
              vectorPerf(it)=results.iterations(it).perf;
              end
              vectorPerfGroup{subjcount}=vectorPerf;
              %% get mean classification accuracy, and SE of that mean
              meanPerf=mean(vectorPerf);
              stePerf=std(vectorPerf)/sqrt(numIterations);
              avgPerfGroup(subjcount)=meanPerf;
              avgPerfSEGroup(subjcount)=stePerf;
              %% get conf matrix
              confmat=multiple_iterations_confusion_amy(results);
              confmatGroup(:,:,subjcount)=confmat;
           
              clearvars results numIterations vectorPerf meanPerf stePerf confmat  
           end 
                end
                
               end
cd([groupmvpadir 'specifics'])      

numsubjectswithROI=sum(hasroi(roin,:))
if numsubjectswithROI>0
%% make summary plot of num voxels in each fold
% % % % scount=0;
% % % % for s=1:numSubj
% % % %     if hasroi(roin,s)
% % % %         scount=scount+1;
% % % %     subplot(pdim,pdim,s);bar(voxelsGroup{scount});ylabel('# voxels');xlabel(subjectlist{s});set(gca,'XTick',[])
% % % %     end
% % % % end
% % % % p=gcf;
% % % % saveas(p, [currentdesc '_voxelcounts.fig']);
% % % % save([currentdesc '_voxelcounts.mat'], 'voxelsGroup');
clear gcf
close all
%% make classification accuracy bar graph (with standard error across folds)
 width=24;
 height=12;
 h1 = figure(); 
set(h1, 'units','inches')
set(h1, 'Position', [10 10 width height])
scount=0;
for s=1:numSubj
    if hasroi(roin,s)
        scount=scount+1;
    subplot(pdim,pdim,s);bar(vectorPerfGroup{scount});ylim([0 1]);ylabel('class perf');xlabel('perf in each fold');title(subjectlist{s});set(gca,'XTick',[])
    refline(0,chance);
    end
end
p=gcf;
saveas(p, [currentdesc '_classPerfInEachFold.fig']);
clear gcf
close all
save([currentdesc '_groupPerfVect.mat'], 'vectorPerfGroup')
%% make classification accuracy subplots (each containing all folds for each subject)
barwitherr(avgPerfSEGroup, avgPerfGroup);ylim([0 1]);ylabel('class perf'); set(gca, 'XTickLabel', subjlabels{roin});
refline(0,chance);
p=gcf;
saveas(p, [currentdesc '_classPerfAcrossFolds.fig']);
clear gcf
close all
save([currentdesc '_groupavgs.mat'], 'avgPerfGroup', 'avgPerfSEGroup')
%% make confusion matrix
scount=0;
 width=24;
 height=12;
 h1 = figure();
set(h1, 'units','inches')
set(h1, 'Position', [10 10 width height])
for s=1:numSubj
    if hasroi(roin,s)
        scount=scount+1;
    subplot(pdim,pdim,s);imagesc(confmatGroup(:,:,scount), [0 1]);ylabel('conf matrix');xlabel(condstring);title(subjectlist{s});set(gca,'XTick',[])
    colormap(hot)
    colorbar
    end
end
 saveas(p, [currentdesc '_confusions.fig']);
clear gcf
close all
save([currentdesc '_confusions.mat'], 'confmatGroup') 

%% make summaries across ROIs
confmatALLROIS(:,:,roicount)=mean(confmatGroup,3);
AvgPerfAllROIS{roicount}=avgPerfGroup;
AvgPerfSEAllROIS{roicount}=avgPerfSEGroup;

clearvars confmatGroup avgPerfGroup avgPerfSEGroup vectorPerfGroup voxelsGroup
cd .. 
cd ..
cd ..
end
            end
foundROIs=roicount;
cd(groupmvpadir)

g=[]
for loc=1:length(localizertasks)
if isempty(g)    
g=findstr(currentdesc, localizertasks{loc});
index=g+length(localizertasks{loc})+1;
end
end

newdesc=currentdesc(index:end);
 width=48;
 height=24;
 h1 = figure(); 
set(h1, 'units','inches')
set(h1, 'Position', [10 10 width height])
for n=1:foundROIs

    roiarray=AvgPerfAllROIS{n};
    if length(roiarray)>1
    [h,p]=ttest(roiarray,chance,0.05,'right'); %% do a ttest comparing mean accuracy across subjects to chance
    else
    h=0;
    end
    if h & length(roiarray)>minsubj2analyze
        barcol='g'; % green if significant
    else
        barcol='r';
    end
subplot(ceil(sqrt(foundROIs)),ceil(sqrt(foundROIs)),n);
barwitherr(AvgPerfSEAllROIS{n}, AvgPerfAllROIS{n}, barcol);
title(roinames{n}); set(gca,'XTick',[]); %set(gca, 'XTickLabel', subjlabels{n});
ylim([0 1]);
if n==1 || mod(n,ceil(sqrt(foundROIs)))==0
ylabel('class perf'); 
end
refline(0,chance)
end
p=gcf;
saveas(p, [newdesc '_classAllROIs.fig']);
clear gcf
close all


%%%%%%%
 width=48;
 height=24;
 h1 = figure(); 
set(h1, 'units','inches')
set(h1, 'Position', [10 10 width height])
kstat=zeros(foundROIs,1);
reject=zeros(foundROIs,1);
pvalue=zeros(foundROIs,1);
for n=1:foundROIs

    roiarray=AvgPerfAllROIS{n};
    if length(roiarray)>4
    [reject(n),pvalue(n),kstat(n)]=lillietest(roiarray); %% do a lilliefors against the null hypothesis that array comes from a normal distribution
    end

subplot(ceil(sqrt(foundROIs)),ceil(sqrt(foundROIs)),n);normplot_AES(roiarray);title(roinames{n}); %set(gca, 'XTickLabel', subjlabels{n});
%ylim([0 1]);
end
p=gcf;
saveas(p, [newdesc '_normplot.fig']);
save([newdesc '_normality_stats.mat'], 'kstat', 'pvalue', 'reject');
clear gcf
close all


 width=48;
 height=24;
 h1 = figure(); 
set(h1, 'units','inches')
set(h1, 'Position', [10 10 width height])
for n=1:foundROIs
    subplot(ceil(sqrt(foundROIs)), ceil(sqrt(foundROIs)),n);imagesc(confmatALLROIS(:,:,n), [0 1]);xlabel(condstring);title(roinames{n});set(gca,'XTick',[]),set(gca,'YTick',[])
    colormap(hot)
    %if mod(n,ceil(sqrt(foundROIs)))==0
    colorbar
    %end
end
p=gcf;
saveas(p, [newdesc '_confusions.fig']);
clear gcf
close all

clearvars confmatALLROIS AvgPerfAllROIS AvgPerfSEAllROIS
            end
        end
end
