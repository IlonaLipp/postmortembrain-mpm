function create_b1_map_with_hmri_toolbox(local_default_file, B1_files,B0_files,parameter_file,outputdir)
    %% makes SPM batch to create SESTE b1+ map with hmri toolbox. previously did have to run whole pipeline, now it just runs the b1 map creation part
    %%% coregistration to PD is done by calculating an anatomical image
    %%% frmo the B1 files, and doing some segmentation (assuming whole
    %%% brain)
    %%% needs parameter file with the correct acquisition parameters
    spm_jobman('initcfg');
    spm('defaults', 'fmri');
    %%% config toolbox
    matlabbatch{1}.spm.tools.hmri.hmri_config.hmri_setdef.customised = {local_default_file};
    %%%
    
%     matlabbatch{2}.spm.tools.hmri.create_mpm.subj.output.outdir = {outputdir};
%     matlabbatch{2}.spm.tools.hmri.create_mpm.subj.sensitivity.RF_none = '-';
%     matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.i3D_EPI.b1input = B1_files;
%     matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.i3D_EPI.b0input = B0_files;
%     matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.i3D_EPI.b1parameters.b1defaults = {parameter_file};
%     matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.MT = {''};
%     matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.PD = PD_files;
%     matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.T1 = {''};
%     matlabbatch{2}.spm.tools.hmri.create_mpm.subj.popup = false;
    
    matlabbatch{2}.spm.tools.hmri.create_B1.subj.output.outdir = {outputdir};
    matlabbatch{2}.spm.tools.hmri.create_B1.subj.b1_type.i3D_EPI.b1input = B1_files;
    matlabbatch{2}.spm.tools.hmri.create_B1.subj.b1_type.i3D_EPI.b0input = B0_files;
    matlabbatch{2}.spm.tools.hmri.create_B1.subj.b1_type.i3D_EPI.b1parameters.b1defaults = {parameter_file};
    matlabbatch{2}.spm.tools.hmri.create_B1.subj.popup = false;
%     
    spm_jobman('run',matlabbatch);

end