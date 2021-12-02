function save_and_brain_mask_files_from_MAPS_structure(MAPS,MPMvol,brain_mask_file,mpm_output_dir, processed_file_string, mpm_run, resolution)
    all_maps = fieldnames(MAPS);
    MPMvol.dt(1) = spm_type('float32');
    for m = 1:length(all_maps)
        mapname = all_maps{m};
        MPMvol.fname = [mpm_output_dir, '/', mapname,'_', processed_file_string, '_', resolution,'_run', sprintf('%.02d',mpm_run),'.nii'];
        spm_write_vol(MPMvol, MAPS.(mapname));
        brain_masked_filename = [mpm_output_dir, '/', mapname,'_', processed_file_string, '_', resolution,'_run', sprintf('%.02d',mpm_run), '_brain_masked.nii'];
        system(['fslmaths ', MPMvol.fname, ' -mas ', brain_mask_file, ' ', brain_masked_filename]);
    end
end