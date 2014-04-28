cd(directory)
p=imagesc([subjectsave.regressors{1}.mat' subjectsave.regressors{2}.mat' subjectsave.regressors{3}.mat']);
    plotptbselectors=gcf;
    saveas(plotptbselectors, [mvpadir printregressor '_srch_' xvalselector '_regressors_noavg.fig']);
    clear gcf
    close(gcf)
    
    p=imagesc([subjectsave.regressors{4}.mat']);
    plotavgselectors=gcf;
    saveas(plotavgselectors, [mvpadir printregressor '_srch_' xvalselector '_regressors_avg.fig']);
    clear gcf
    close(gcf)
    
    notaveraged=[];
    averaged=[];
    for x=1:length(subjectsave.selectors)
       name=subjectsave.selectors{x}.name;
       if isempty(strfind(name, 'avg'))
           notaveraged=[notaveraged (subjectsave.selectors{x}.mat/max(subjectsave.selectors{x}.mat))'];
       else 
           averaged=[averaged (subjectsave.selectors{x}.mat/max(subjectsave.selectors{x}.mat))'];
       end
    end
    
    p=imagesc(notaveraged);
    noavgselectors=gcf;
    saveas(noavgselectors, [mvpadir printregressor '_srch_' xvalselector 'selectors_noavg.fig']);
    clear gcf
    close(gcf)
    
    p=imagesc(averaged);
    avgselectors=gcf;
    saveas(avgselectors, [mvpadir printregressor '_srch_' xvalselector 'selectors_avg.fig']);
    clear gcf
    close(gcf)
    
     
    %turns out you'll use a lot of space even if you just save the full
    %subjstructure for every mask. copy and delete out masks from structure
    %that is saved to dir and just save it once
    subjectsave2dir=subjectsave;
    subjectsave2dir.masks=[]; % to replace a given mask m into a structure later on just set subjectsave2dir.masks=savemasks{m}
    
    %% save the important stuff
    %save([mvpadir mask_description '_' printregressor '_' xvalselector '_classification.mat'], 'results');
    save([mvpadir printregressor '_' xvalselector '_srch_subjstructure.mat'], 'subjectsave2dir');
