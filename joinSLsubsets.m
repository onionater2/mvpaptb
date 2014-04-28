function joinSLsubsets(study, resultsdir, subjectlist)
rootdir=['/mindhive/saxelab2/' study '/']
deleteimgs='YES'; %'Yes; will delete images, 'No' will keep them (unless you set it to ask you at each mask, which will override this)
askbeforedelete=0; %will ask after creating each union whether you want to delete subset images
for s=1:length(subjectlist)
    subject=subjectlist{s}
    subjSLdir=[rootdir, subject '/mvpa_ptb/' resultsdir];
    cd(subjSLdir)
    imgfiles=dir('*.img');
    for f=1:length(imgfiles)
        filename=imgfiles(f);
        filename=filename.name;
        imgfilescell{f}=filename(1:end-7);
    end
    imgtypes=unique(imgfilescell);

    for i=1:length(imgtypes)
        imagetype=imgtypes{i}
        subsetcheck=imagetype(end-5:end);
        theseimages=[];
        data=[];
        uniondata=[];
        if strcmp('subset', subsetcheck)
        theseimages=dir([imagetype '*.img']);
        for img=1:length(theseimages)
        subsetimg=theseimages(img).name;
        file=spm_vol(subsetimg);
        data(:,:,:,img)=spm_read_vols(file);
        
        %delete this img
        end
        uniondata=sum(data,4);
        writeTemplate=file(1); % getting template from a single fold classification .img, keeping pinfo
        classificationpinfo=writeTemplate.pinfo
        writeTemplate.dt = [spm_type('float64') spm_platform('bigend')];
        writeTemplate.fname = [imagetype(1:end-6) 'union.img'];
        unionoutput.fname =[imagetype(1:end-6) 'union.img'];
        unionoutput=spm_create_vol(writeTemplate);
        unionoutput = spm_write_vol(unionoutput, uniondata);        
        
        if askbeforedelete
        deleteimgs = questdlg([imagetype ' union image created.  ','Delete subset images?']);
        end
        if strcmpi(deleteimgs,'Yes');
        for img=1:length(theseimages)
        subsetimg=theseimages(img).name
        subsethdr=[subsetimg(1:end-3) 'hdr'];
        delete(subsetimg)
        delete(subsethdr)
        end
        end

    end
    end
end
end
    