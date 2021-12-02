function file = get_first_nifti_file_from_folder(folder)
    all_files = dir([folder,'/*.nii']);
    file = [all_files(1).folder,'/',all_files(1).name];
end