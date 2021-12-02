function process_mp2rage(param, inv2_file, uni_file, t1_file, outputdir, bm_thr_mp2rage, run)

    basename = [outputdir,'/MP2RAGE_repeat_',sprintf('%.02d',run),'_'];

    signal_mask_file = [basename,'signal_mask.nii'];

   system(['fslmaths ',inv2_file,' -thr ',num2str(bm_thr_mp2rage),' -bin ',signal_mask_file]); %%% threshold based on signal intensity of raw file

   if exist(t1_file) ~= 0%%% not there for all scans   
     
      t1_file_sm = [basename, 'T1_calculated_from_scanner_signal_masked.nii'];
      system(['fslmaths ', t1_file,' -mul ',signal_mask_file, ' ',t1_file_sm]);
      
      r1_file_sm = [basename, 'R1_calculated_from_scanner_signal_masked.nii'];
      system(['fslmaths ',signal_mask_file,' -div ',t1_file_sm, ' -mul 1000 ',r1_file_sm]);
      
   end
   if exist(param.qmrlabpath) ~= 0
        Guiopts = load(param.Guiopts);
        FitResults = calculate_t1_with_qmri_lab(param.qmrlabpath, uni_file, Guiopts);
        spmvol = spm_vol(uni_file);
        spmvol.dt(1) = spm_type('float32');
        spmvol.fname = [basename,'T1_calculated_with_qmrlab_signal_masked.nii'];
        spm_write_vol(spmvol, FitResults.T1 * 1000); %%% in MS
%         if strcmp(param.useforT1ref, 'qmrlab')
%            system(['fslstats ', spmvol.fname, ' -P 50 >> ', spmvol.fname(1:end-4),'.txt']);
%         end
        system(['fslmaths ', spmvol.fname,' -mul ',signal_mask_file, ' ',spmvol.fname]);
        spmvol.fname = [basename,'R1_calculated_with_qmrlab_signal_masked.nii'];
        spm_write_vol(spmvol, FitResults.R1);
        system(['fslmaths ', spmvol.fname,' -mul ',signal_mask_file, ' ',spmvol.fname]);

   end
   
end