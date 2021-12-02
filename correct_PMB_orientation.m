function correct_PMB_orientation(input, output, permutation_order, flipdir)

    %%% permutes axis to make orientation labels in FSL be correct (this is
    %%% needed because postmortem brains do not lie in the scanner in the
    %%% same orientation as in vivo people do. 
    %%% needs permutation order (e.g. [3 2 1] to swap x and z axis) and if
    %%% you wanted to flip any axis (for whatever reason) you can indicate
    %%% the axes to flip as well. input is nifti filename, output too.
    
    test = spm_vol(input);
    img = spm_read_vols(test);

    %test.mat = test.mat([permutation_order 4],:); %%% we actually do not
    %need to correct sform in image

    test.dim = test.dim([permutation_order]);
    imgnew = permute(img,[permutation_order]);
    
    %%% flip dimensions if needed
    if length(flipdir) > 0
        for toflip = flipdir
            imgnew = flipdim(imgnew,[toflip]);
        end
    end
    
    test.fname = output;

    spm_write_vol(test,imgnew);

end

% for emma 190326 permutation order needs to be [3 2 1], for groat 191031
% is also different
% input='/data/pt_02101/preprocessed/006_C_W_NORTHEAST1_TAI_NE/mr/191019_Magnetom_7T_32Ch_WB/MPMs/R1_0p4_run01_brain_masked.nii'
% output='/home/raid2/lippi/Documents/scripts/Primate_resources/SamVickery_Templates/amended/Northeast_191019_R1_0p4_run01_brain_masked_orientation_corrected.nii'
% permutation_order = [1 3 2];

% permutation_order = [3 2 1]; %%% don't ask why
% permutation_order = [2 1 3]
% permutation_order = [3 1 2]
