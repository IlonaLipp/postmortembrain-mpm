%file_B1 = '/data/pt_02101/preprocessed_data/ALBI_190314_7T_32Ch/MPMs/B1map_highres_0p4_run3.nii';
%file_BS_B1 = '/data/pt_02101/preprocessed_data/ALBI_190314_7T_32Ch/MPMs/BS_B1map_highres_0p4_run3.nii';

function compare_b1_maps_with_BA(mapfile1, mapfile2, figurefilename)
    %%% this function compares the values of two B1+ maps, whereby the mapfiles need to be coregistered
    %%% it outputs a blant-altman plot

    pathname = fileparts(mapfile1)
    
    B1map1 = spm_vol(mapfile1); 
    B1map1 = spm_read_vols(B1map1);
    B1map1_vec = B1map1(:);

    B1map2 = spm_vol(mapfile2); 
    B1map2 = spm_read_vols(B1map2);
    B1map2_vec = B1map2(:);
    
    mean_vals = mean([B1map1_vec,B1map2_vec],2);
    %idx = find(mean_vals > 50 & mean_vals < 110); %%% can change range of considered values to make plot easier to interpret
    idx = find(mean_vals > 0);
    
    mean_vals = mean([B1map1_vec(idx),B1map2_vec(idx)],2);
    diff_vals = B1map1_vec(idx) - B1map2_vec(idx);
    mean(mean_vals)
    mean(diff_vals)
    %[N,C] = hist3([mean_vals,diff_vals],[100,50]);
    
    %[N,C] = hist3([mean_vals,diff_vals],{linspace(50,125,100),linspace(-15,15,100)});
    [N,C] = hist3([diff_vals,mean_vals],{linspace(-15,15,100),linspace(50,125,100)});
    %toplot = [N(:,end:-1:1)];
    toplot = [N(end:-1:1,:)];
    imagesc(toplot)
    xticks(1:4:length(C{2}))
    xticklabels(round(C{2}(1:4:end),2))
    xtickangle(90)
    set(gca,'FontSize',11)
    yticks(1:10:length(C{1}))
    yticklabels(round(C{1}(end:-10:1),2))

end