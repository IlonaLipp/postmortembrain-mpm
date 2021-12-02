function series_numbers = dir_names_to_series_numbers(dirnames)
    %%% function to read in output of matlabs dir function over a scan
    %%% file dirctory and returns a vector of series numbers
    if length(dirnames) > 0
        for d = 1:length(dirnames)
           curr_name = dirnames(d).name;
           split = strsplit(curr_name,'S');
           split2 = strsplit(split{2},'_');
           series_numbers(d) = str2double(split2{1});
        end
    else
        series_numbers = [];
    end
    series_numbers = sort(series_numbers);
end