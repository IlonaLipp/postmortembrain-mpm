function MPMcalc_basic_with_hmri_toolbox_b1corr(hmri_processing_parameters,outdir, MTwimages, PDwimages, T1wimages, b1map)

    %%% make b1map structure
    b1maps{1} = b1map;
    b1maps{2} = b1map;
    
    spm_jobman('initcfg');
    spm('defaults', 'fmri');
    matlabbatch{1}.spm.tools.hmri.hmri_config.hmri_setdef.customised = {hmri_processing_parameters};
    matlabbatch{2}.spm.tools.hmri.create_mpm.subj.output.outdir = {outdir};
    matlabbatch{2}.spm.tools.hmri.create_mpm.subj.sensitivity.RF_none = '-';
    matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.pre_processed_B1.b1input = b1maps';
    matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.pre_processed_B1.scafac = 1;
    matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.MT = MTwimages;
    matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.PD = PDwimages;
    matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.T1 = T1wimages;
    matlabbatch{2}.spm.tools.hmri.create_mpm.subj.popup = false;
    spm_jobman('run',matlabbatch);
    
end