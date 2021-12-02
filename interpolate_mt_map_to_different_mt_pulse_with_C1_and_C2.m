function newmap = interpolate_mt_map_to_different_mt_pulse_with_C1_and_C2(mtmap, C1, C2, nominal_fa_rad, desired_fa_rad)
    %%% make fake b1+ map
    fake_b1 = 100 * nominal_fa_rad / desired_fa_rad;
    %fake_b1 = 100 * desired_fa_rad / nominal_fa_rad; %%% trying if this is better
    fake_b1map = fake_b1 * ones(size(mtmap));
    %newmap = b1_correct_mt_map(mtmap, fake_b1map, C); %%% -C because relationship is in the opposite direction now
    map_of_ones = ones(size(mtmap));
    newmap = mtmap ./ ((map_of_ones - C1 + C2) ./ (map_of_ones - C1 .* fake_b1map/100 + C2 .* (fake_b1map/100).^2  ));
    
end