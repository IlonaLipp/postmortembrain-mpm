function create_afi_b1_map_with_hmri_toolbox(local_default_file,param,B1_files,PD_files,outputdir)
    %%% this function creates an SPM batch to calculate AFI b1+ maps, with
    %%% the hMRI toolbox

    spm_jobman('initcfg');
    spm('defaults', 'PET');
    addpath(genpath([param.bashscriptfolder,'/hMRI_toolbox_patch/'])); %%% patch that i created that also saves some intermediate steps of the AFI map processing not just the end step so that i can use ...
    %%% my own processing pipeline on them
    matlabbatch{1}.spm.tools.hmri.hmri_config.hmri_setdef.customised = {local_default_file};

%     matlabbatch{2}.spm.tools.hmri.create_mpm.subj.output.outdir = {outputdir};
%     matlabbatch{2}.spm.tools.hmri.create_mpm.subj.sensitivity.RF_us = '-';
%     matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.i3D_AFI.b1input = B1_files;
%     matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.i3D_AFI.b1parameters.b1metadata = 'yes';
%     matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.MT = '';
%     matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.PD = PD_files;
%     matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.T1 = '';
%     matlabbatch{2}.spm.tools.hmri.create_mpm.subj.popup = false;
    
    matlabbatch{2}.spm.tools.hmri.create_B1.subj.output.outdir = {outputdir};
    matlabbatch{2}.spm.tools.hmri.create_B1.subj.b1_type.i3D_AFI.b1input = B1_files;
    matlabbatch{2}.spm.tools.hmri.create_B1.subj.b1_type.i3D_AFI.b1parameters.b1defaults = {local_default_file};
    matlabbatch{2}.spm.tools.hmri.create_B1.subj.popup = false;
    spm_jobman('run',matlabbatch);
    
    
    
end