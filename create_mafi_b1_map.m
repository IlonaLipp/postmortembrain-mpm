function create_mafi_b1_map(file1, file2, desired_fa, tr_diff, outputfile)
    S1 = spm_read_vols(spm_vol(file1));
    S2 = spm_read_vols(spm_vol(file2));            
    %B1rel = Saskia_AFI_B1Map_SH(S1,S2,tr_diff,desired_fa);
    %%% copied from Saskia_AFI_B1Map_SH
    r = S1./S2;
    n = tr_diff;
    nominalAlpha = deg2rad(desired_fa);
    B1rel = acos((r*n-1)./(n-r))./nominalAlpha;
    %%% write out
    file_templ = spm_vol(file1);
    file_templ.fname = outputfile;
    file_templ.dt(1) = spm_type('float32');
    spm_write_vol(file_templ, B1rel * 100);
end