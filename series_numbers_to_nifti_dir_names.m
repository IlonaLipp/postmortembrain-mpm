function foldernames = series_numbers_to_nifti_dir_names(parent_dir, series_numbers)
    %%% this function takes a vector of series numbers and creates a list
    %%% of directory names (from the nfiti folder) from them
    if length(series_numbers) > 0
        for d = 1:length(series_numbers)
            diroutput = dir([parent_dir,'/*',sprintf('%04d',series_numbers(d))]);
            foldernames(d).folder = parent_dir; 
            foldernames(d).name = diroutput(1).name; 
        end
    else
        foldernames = [];
    end
end