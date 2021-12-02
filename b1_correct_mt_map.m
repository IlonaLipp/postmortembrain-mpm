function newmap = b1_correct_mt_map(mtmap, b1map, C)
%%% C = calibration constant C from Gunther's approach
%%% mt map is map to be calibrated
%%% b1 map in %
%%% newmap = b1 corrected mtmap
    %b1_to_corr = 100*ones(size(mtmap)) - b1map; 
    %%% newmap = mtmap + mtmap .* b1_to_corr .* slope; %%% multiplies the unnormalised slope with b1+
    map_of_ones = ones(size(mtmap));
    newmap = mtmap .* ((map_of_ones - C) ./ (map_of_ones - C .* b1map/100));
end