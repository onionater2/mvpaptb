function [dataall, datalabel, subjectname]= compare_mvpastreams(folder, subjectlist)
%created by AES to compare results of different mvpa analysis streams

rootdir='/mindhive/saxelab2/EIB/'
fo=[rootdir, 'mvpaptb/compare_streams_' datestr(date) '.csv'];
f=fopen(fo,'w');

hardcodedrois= {
    'rinsula_wfu_xyz_group',...
    'rvSTR_reward_xyz_group',...
    'vmPFC_reward_xyz_group',...
    'right_ant_temporal_xyz_group',...
    'ramygdala_wfu_xyz_group',...
    'MPFC_peelenpeak_xyz_group',...
    'lSTS_peelenpeak_xyz_group',...
    'rSTS_peelenflip_xyz_group',...
    'lvSTR_reward_xyz_group',...
    'linsula_wfu_xyz_group',...
    'left_ant_temporal_xyz_group',...
    'lamygdala_wfu_xyz_group',...
    'RTPJ_tomloc_ind',...
    'LTPJ_tomloc_ind',...
    'RSTS_tomloc_ind',...
    'LSTS_tomloc_ind',...
    'PC_tomloc_ind',...
    'MMPFC_tomloc_ind',...
    'DMPFC_tomloc_ind',...
    'rpSTS_BDbiomot_ind',...
    'rFFA_kanparcelFaceObj_EmoBioLoc',...
    'lFFA_kanparcelFaceObj_EmoBioLoc',...
    'rSTS_kanparcelFaceObj_EmoBioLoc',...
    'lSTS_kanparcelFaceObj_EmoBioLoc',...
    'rOFA_kanparcelFaceObj_EmoBioLoc',...
    'lOFA_kanparcelFaceObj_EmoBioLoc',...
    'rLOC_foundObjFace_EmoBioLoc',...
    'lLOC_foundObjFace_EmoBioLoc'};

numROIs= length(hardcodedrois);

for s=1:length(subjectlist)
    subject=subjectlist{s}
mvpadir=[rootdir subject '/' folder '/']
cd(mvpadir)


disc=load('discriminations.mat');
names=fieldnames(disc);
numdisc=length(names)


contents=dir;
if s==1
numstuff=length(contents);
countstreams=0;
for t=3:numstuff 
    if isdir(contents(t).name) countstreams=countstreams+1; end
end
dataall(:, :,length(subjectlist))=zeros(countstreams, 2*numdisc*numROIs);
datalabel=cell(1,2*numdisc*numROIs);
end

count=0;

%
for dindex=1:numdisc
    for roiindex=1:numROIs
        for selindex=1:2
   dataindex=((dindex-1)+(dindex+(selindex-1))+(roiindex-1)*(numdisc*2));
   sels{1}='runsel';
   sels{2}='othersel';
   thissel=sels{selindex};
   datalabel{dataindex}=[hardcodedrois{roiindex} '_' names{dindex} '_' thissel];
        end
    end
end

for x=3:length(contents)
   
    if isdir(contents(x).name)
        cd(contents(x).name);
        count=count+1;

        subjectname(s).analysisname{count}=contents(x).name
        
        classfiles=dir('*classification.mat');
        
        for c=1:length(classfiles)
        strpat=classfiles(c).name(1:12);
        
        for d=1:numdisc
            finddisc = strfind(classfiles(c).name, names{d});
            if ~isempty(finddisc)
                dindex=d;
                discname=names{d};
                findsel=strfind(classfiles(c).name, 'run');
                if ~isempty(findsel)
                    discaddindex=0;
                    sel='runsel';
                else
                    discaddindex=1;
                    sel='othersel';
                end
                
            end
        end
        
        for r=1:numROIs
            roi=hardcodedrois{r};
           if strcmp(roi(1:12), strpat)
               roi;
               roiindex=r;
           end
        end
        
        dataindex=((dindex-1)+(dindex+discaddindex)+(roiindex-1)*(numdisc*2));
        load(classfiles(c).name);
        numfolds=length(results.iterations);
        if numfolds>2
            tempdata=results.total_perf;
        else if numfolds==2
        for n=1:numfolds
            tempiterdata(n)=results.iterations(n).perf;
            tempdata=mean(tempiterdata);
        %load some relevant data
        end
            end
        end
        dataall(count,dataindex,s)=tempdata;
        end
        cd ..
        
        % print it out
        
        name=subjectname(s).analysisname{count};
        datavect=dataall(count,:,s);
        numberCol=size(datavect,2);
        inputstring=['%d'];
        inputs=s;
        for d=1:numberCol
        inputstring=[inputstring ' %d'];
        inputs=[inputs; datavect(1,d)];
        end
        fprintf(f, inputstring, [inputs']);
        fprintf(f,'\n');
        
        
    
    
end

end
end
fclose(f);
end