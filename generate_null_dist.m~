function generate_null_dist(study, subjectlist, mvparesults, task, runsincluded, numBootstraps, desiredSig, observedacc)

numSubj=length(subjectlist);
pdim=ceil(sqrt(numSubj));
indMVPAdir='mvpa_ptb/';
localizertasks={'tomloc', 'EmoBioLoc'}

firstsubj=subjectlist{1}
lastsubj=subjectlist{numSubj}
subjectrange=['_subj' firstsubj(end-1:end) 'to' lastsubj(end-1:end)];

rootdir=['/mindhive/saxelab2/' study '/'];
cd(rootdir)
mkdir([rootdir 'mvpaptb/' mvparesults, subjectrange '/']);
subjmvpa=[indMVPAdir mvparesults '/'];
groupmvpadir=[rootdir 'mvpaptb/' mvparesults, subjectrange '/'];


%roilist={'whole_brain_skullstripped'};
roilist={
    'RTPJ_tomloc',...
    'LTPJ_tomloc',...
    'RSTS_tomloc',...
    'LSTS_tomloc',...
    'PC_tomloc',...
    'MMPFC_tomloc',...
    'DMPFC_tomloc',...
    'MPFC_peelenpeak_xyz',...
    'lSTS_peelenpeak_xyz',...
    'rpSTS_BDbiomot',...
    'rFFA_kanparcel_EmoBioLoc',...
    'lFFA_kanparcel_EmoBioLoc',...
    'rSTS_kanparcel_EmoBioLoc',...
    'lSTS_kanparcel_EmoBioLoc',...
    'rOFA_kanparcel_EmoBioLoc',...
    'lOFA_kanparcel_EmoBioLoc',...
    'rLOC_found_EmoBioLoc',...
    'lLOC_found_EmoBioLoc',...
    'rinsula_wfu_xyz',...
    'rvSTR_reward_xyz',...
    'vmPFC_reward_xyz',...
    'right_ant_temporal_xyz',...
    'ramygdala_wfu_xyz',...
    'lvSTR_reward_xyz',... 
    'linsula_wfu_xyz',...
    'left_ant_temporal_xyz',...
    'lamygdala_wfu_xyz',...
    'rSTS_peelenflip_xyz'};



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
            foundtemplate=0;
            for fs=1:numSubj
            for fr=1:length(roilist)       
            lookin=[rootdir subjectlist{fs} '/' subjmvpa roilist{fr} '*' desc '*classification_iter1.mat']; 
            templateDES=dir(lookin); % to identify num of disc variants for the descrimination
            if ~isempty(templateDES)
                foundtemplate=1;
                break % from roi loop
            end
            end
            if foundtemplate
                break %from subj loop
            end
            end
            numDES=size(templateDES);
            for dse=1:numDES(1)
            tempdesc=templateDES(dse).name;
            roicount=0;
            hasroi=zeros(numROIs, length(subjectlist));
            for roin=1:numROIs
               roiname=roilist{roin}; 
               index=strfind(tempdesc, desc);
               tempdesc=tempdesc(index:end);
               disp(['working on discrimination ' tempdesc ' in roi ' roilist{roin}]);
               subjindices=zeros(numSubj,1);
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
                  
                    subjdescs=dir([resultsdir roiname '*' tempdesc(1:end-5) '*']);
                    if isempty(subjdescs)
                        disp(['subject ' subject ' doesnt have that discrimination for roi ' roiname])
                        hasroi(roin,subj)=0;
                    else
                    permiterations=length(subjdescs);
                    load(subjdescs(1).name);
                    numIterations=length(results.iterations);
                    if ~exist('avgPerfGroup') avgPerfGroup=zeros(permiterations, numSubj); end
                    if ~exist('vectorPerfGroup') vectorPerfGroup=zeros(permiterations, numIterations, numSubj); end
                    if ~exist('avgPerfSEgroup') avgPerfSEGroup=zeros(permiterations, numSubj); end 

                    for rep=1:permiterations 
                    subjdesc=subjdescs(rep);
                    currentdesc=subjdesc.name;
                    currentdesc=currentdesc(1:end-19);
                        if hasroi(roin,subj)==0
                        hasroi(roin,subj)=1;
                        end
                        subjindices(subj)=1;
                

              %% get vector of classification accuracies across folds
              newit=0;
              newit=load(subjdesc.name);
              numConds=length(unique(newit.results.iterations(1).perfmet.desireds));
              chance=1/numConds;
              for it=1:numIterations
              vectorPerf(rep,it)=newit.results.iterations(it).perf;
              end
              vectorPerfGroup(rep,:, subj)=vectorPerf(rep,:);
              %% get mean classification accuracy, and SE of that mean
              meanPerf(rep,1)=mean(vectorPerf(rep,:));
              stePerf(rep,1)=std(vectorPerf(rep,:))/sqrt(numIterations);
              avgPerfGroup(rep,subj)=meanPerf(rep);
              avgPerfSEGroup(rep,subj)=stePerf(rep);
              %% get conf matrix

              end
              save([currentdesc '_classification_distribution.mat'], 'meanPerf', 'stePerf', 'vectorPerf')
              clearvars results numIterations vectorPerf meanPerf stePerf   
           end 
                end               
               end
      

numsubjectswithROI=sum(hasroi(roin,:))
if numsubjectswithROI>0


%%bootstrap group-level accuracies from these distributions
insubjects=find(subjindices);
for b=1:numBootstraps
for i=1:length(insubjects)
    thisubjvector=avgPerfGroup(:,insubjects(i));
    randomacc=randi([1,permiterations]);
    groupvector(i)=thisubjvector(randomacc);
end
groupacc(b)=mean(groupvector);
end

sortedgroupacc=sort(groupacc);
numaccs=length(groupacc);
countfalsepositives=0;
maxFPelements=ceil(numaccs*desiredSig);
minAcc=sortedgroupacc(end-maxFPelements);
for x=1:numaccs
   if groupacc(x)>observedacc 
       countfalsepositives=countfalsepositives+1;
   end
end
step=.005;
meanacc=mean(groupacc);

thisacc=groupacc*0; thisacc(1)=observedacc;
thisacchist=histc(thisacc, 0:step:1);
thisacchist(1)=0; %(don't count the zeros)
acchist=histc(groupacc, 0:step:1);
countFPs=0;
for bin=length(acchist):-1:1
    countFPs=countFPs+acchist(bin);
    if countFPs>=maxFPelements
        threshold=bin;
        righttail=bin+1;
        break
    end
end
dimensions=find(acchist);
if acchist(find(thisacchist))==0
acchist(find(thisacchist))=1;
end
h=bar(diag(acchist),1, 'stacked');
set(gca, 'xlim', [min(dimensions)-ceil(.05/step) max(dimensions)+ceil(.05/step)])
bindex=find(thisacchist);
minaxis=min(dimensions)-ceil(.05/step);
maxaxis=max(dimensions)+ceil(.05/step);
minaxisval=minaxis*step;
maxaxisval=maxaxis*step;
axisrange=maxaxis-minaxis;
axislabelrange=maxaxisval-minaxisval;
axisstep=axisrange/6;
axislabelstep=axislabelrange/6;
set(gca, 'xtick', [minaxis:axisstep:maxaxis])
axisticksnums=[minaxisval:axislabelstep:maxaxisval];
set(gca, 'xticklabels', axisticksnums)
set(h(minaxis:righttail),'facecolor','b');
set(h(threshold),'facecolor',[.3 .7 .6]);
set(h(righttail:maxaxis),'facecolor','g');
set(h(bindex),'facecolor','r');

statsoutput.currentdesc=currentdesc
statsoutput.numaccs=numaccs;
statsoutput.meanacc=meanacc;
statsoutput.numfps=countfalsepositives;
statsoutput.pvalue=countfalsepositives/numaccs;
statsoutput.desiredSig=desiredSig;
statsoutput.minaccAtDesiredSig=minAcc
  
%% make group summaries
cd(groupmvpadir)
%saveas(gcf, [currentdesc '_permutation_hist.fig'])
close(gcf)
save([currentdesc '_permutation_analysis.mat'], 'avgPerfGroup', 'avgPerfSEGroup', 'vectorPerfGroup', 'groupacc', 'statsoutput')

    

clearvars avgPerfGroup avgPerfSEGroup vectorPerfGroup thisubjecvector groupvector groupacc maxFPelements minAcc statsoutput
cd ../../.. 
end
            end

            end
        end
end
