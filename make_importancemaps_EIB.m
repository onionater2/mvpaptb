function make_importancemaps_EIB(study, task, mvparesultsfolder, subjectlist, runsincluded)
%written by AES 7/31/13 to create importance maps based on voxel weights.
%calls mvpaptb function interpret_weights.m
rootdir='/mindhive/saxelab2/';
studydir=[rootdir study '/' ];

for s=1:length(subjectlist)
    subject=subjectlist{s}
    subjmvpadir=[studydir subject '/mvpa_ptb/' mvparesultsfolder]
    maskdir=[studydir subject '/3danat/'];
    cd(subjmvpadir)
    descmats=dir('*classification*mat')
    numdesc=length(descmats)
    
    % figure out possible discriminations
    template=load([studydir subject '/mvpa_ptb/discriminations_' task '_' runsincluded]);
    possdescs=fieldnames(template);
    %templatesubjinfo=load([rootdir subjectlist{1} '/' indMVPAdir 'subjinfo_' task '_' runsincluded]);
    
    for d=1:numdesc
        desc=descmats(d).name;
        if strcmp(desc(1:10), 'whole_brai')
        desc=desc(1:end-19)
        for thisd=1:length(possdescs)
        findit=strfind(desc,possdescs{thisd});
        if ~isempty(findit)
            break
        end
        end
        thisdescsubjstruct=[desc(findit:end) '_subjstructure.mat'];
        load(thisdescsubjstruct);
        load(descmats(d).name);
        thissubject=subjectsave2dir(1)
        thissubject = load_analyze_mask(thissubject,'whole_brain_skullstripped',[maskdir 'grey_matter_mask.img']);
        thissubject=interpret_weights_aes(thissubject,results) %existin interpret_weights only works for backprop. intend to make my own but haven't done anything to it yet
    end
    end
end



end