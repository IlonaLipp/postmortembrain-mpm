function cell_structure = convert_dir_output_to_cell_structure(filenames)
    %%% helper function to create a cell structure of file names from the
    %%% output of matlab's "dir" function
    if length(filenames) > 0
        for file = 1:length(filenames)
            cell_structure{file} = [filenames(file).folder,'/',filenames(file).name];
        end
        cell_structure = cell_structure';
    else
        cell_structure = {};
    end
    