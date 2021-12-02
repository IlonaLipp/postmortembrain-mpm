function map = calculate_receive_bias_map(transmit_nifti_folder, receive_nifti_folder, output_filename)
 
 %%% load in all echos
 rec = make_mpm_file_structure_from_filenames(dir([receive_nifti_folder,'/*nii']));
 tra = make_mpm_file_structure_from_filenames(dir([transmit_nifti_folder,'/*nii']));
 
 for vol = 1:(min([5,length(rec.TEs)]))
    recmap = rec.data(:,:,:,vol);
    rec_vox(:,vol) = recmap(:);
    tramap = tra.data(:,:,:,vol);
    tra_vox(:,vol) = tramap(:);
    %%% check how similar contrast is between transmit echo and first echo
    %%% of receive, assuming receive bias scans all look good
    spatial_similarity = corrcoef(rec_vox(:,1), tra_vox(:,vol));
    all_similarity(vol) = spatial_similarity(1,2);
 end

 %%% for each included echo, calculate a bias map
 included_echos = find(all_similarity > 0.8); %%% completely arbitrary threshold
 for e = 1:length(included_echos)
     echo = included_echos(e);
     biasmap(:,:,:,e) = rec.data(:,:,:,echo) ./ tra.data(:,:,:,echo);
 end
 
 %%% average across all considered echos
 map = median(biasmap,4);
 
 %%% save
 templ_file = dir([receive_nifti_folder,'/*nii']); templ_file = [templ_file(1).folder,'/',templ_file(1).name];
 templ = spm_vol(templ_file);
 templ.fname = output_filename;
 templ.dt(1) = spm_type('float32');
 spm_write_vol(templ, map);
end