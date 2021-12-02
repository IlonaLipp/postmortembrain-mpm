function convert_dicoms_to_nifti_connectom(dicom_data_dir, nifti_data_dir)
    filenames = dir([dicom_data_dir]);
        if length(filenames) > 3
        count = 0;
        for file = 1:length(filenames)
            if length(filenames(file).name) > 10
                count = count + 1;
                ss = strsplit(filenames(file).name,'.00');
                seq(count) = str2double(ss{2}(1:2));
                filename_structure_all{count} = [filenames(file).folder,'/',filenames(file).name];
            end
        end
        all_s = unique(seq);
        for s = 1:length(all_s)
            idx = find(seq == all_s(s));
            nifti_there = dir([nifti_data_dir,'/*',sprintf('%.4d',s),'/*.nii']);
            if length(nifti_there) == 0
                filename_structure = filename_structure_all(idx);
                spm_jobman('initcfg');
                spm('defaults', 'FMRI');
                dicom2nifti_batch{1}.spm.tools.hmri.dicom.data = filename_structure';
                dicom2nifti_batch{1}.spm.tools.hmri.dicom.root = 'series';
                dicom2nifti_batch{1}.spm.tools.hmri.dicom.outdir = {nifti_data_dir};
                dicom2nifti_batch{1}.spm.tools.hmri.dicom.protfilter = '.*';
                dicom2nifti_batch{1}.spm.tools.hmri.dicom.convopts.format = 'nii';
                dicom2nifti_batch{1}.spm.tools.hmri.dicom.convopts.metaopts.mformat = 'sep';
                dicom2nifti_batch{1}.spm.tools.hmri.dicom.convopts.icedims = 0;
                try
                    spm_jobman('run',dicom2nifti_batch);
                end
            end
        end
    end
end