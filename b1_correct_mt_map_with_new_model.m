function newmap = b1_correct_mt_map_with_new_model(mtmap, b1map, C1, C2)
%%% Cs = calibration constant Cs from Gunther's new model approach
%%% mt map is map to be calibrated
%%% b1 map in %
%%% newmap = b1 corrected mtmap
    %b1_to_corr = 100*ones(size(mtmap)) - b1map; 
    %%% newmap = mtmap + mtmap .* b1_to_corr .* slope; %%% multiplies the unnormalised slope with b1+
    map_of_ones = ones(size(mtmap));
    %newmap = mtmap .* ((map_of_ones - C) ./ (map_of_ones - C .* b1map/100));
    newmap = mtmap .* ((1 - C1 + C2) ./ (1 - C1 .* b1map/100 + C2 .* ((b1map/100).^2)));

end