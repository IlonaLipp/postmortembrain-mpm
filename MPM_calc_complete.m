function MPM_calc_complete(param)
    %%% this is the processing pipeline, which will do everything, given
    %%% specifications in the param structure

    %%% specify outputfilenames for QC stuff
    pdfoutputfilename = [param.processed_data_folder,'/',param.unique_id,'_session_parameters.pdf'];
    temperature_traceoutputfile = [param.processed_data_folder,'/',param.unique_id,'_temperature_trace.mat'];
    temperature_figureoutputfile = [param.processed_data_folder,'/',param.unique_id,'_temperature_trace.png'];
    frequencyplot_outputfilename = [param.processed_data_folder,'/',param.unique_id,'_frequency_trace.png'];
    seriesinfo_outputfilename = [param.processed_data_folder,'/',param.unique_id,'_seriesinfo.mat'];
    
    %% dicom conversion
    %%% create a folder for nifti files if it does not exist:
    if exist(param.raw_data_folder_nifti) == 0
        system(['mkdir ', param.raw_data_folder_nifti]);
    end
    %%% run function 'convert_dicoms_to_nifti', which converts dicoms in folders that have empty corresponding nifti folders
    if param.redo_nifti_conversion == 1
        if param.do_nifti_conversion_connectom_style == 1 %%% if the expected data structure is connectom style
            convert_dicoms_to_nifti_connectom(param.raw_data_folder_dicom, param.raw_data_folder_nifti);
        else %%% if the expected data structure is 7T style
            convert_dicoms_to_nifti(param.raw_data_folder_dicom, param.raw_data_folder_nifti, param.do_not_check_niftis);
        end
    end
    
    %% find nifti folders
    %%% takes as an input the sequence numbers in the param structure and
    %%% finds the respective folder name for that series
    all_sequences = fieldnames(param.folder_series_numbers);
    for s = 1:length(all_sequences)
        seqname = all_sequences{s};
        try
            folders.(seqname) =  series_numbers_to_nifti_dir_names(param.raw_data_folder_nifti,param.folder_series_numbers.(seqname));
        end
    end

    %% make output directiory 
    if exist(param.processed_data_folder) == 0
        system(['mkdir ', param.processed_data_folder]);
    end

    %% create scan summary pdf and frequency drift plot and save scan times 
    if exist(pdfoutputfilename) == 0 
        if exist(seriesinfo_outputfilename) == 0
            try
                seriesinfo = analyse_scan_parameters_from_json_files_for_a_scan_session(param,param.raw_data_folder_nifti, pdfoutputfilename, frequencyplot_outputfilename)
                save(seriesinfo_outputfilename,'seriesinfo');
            end
        end
    end

    %% try to make temperature plot (if a temperature file is provided in param)
    if exist(temperature_traceoutputfile) == 0
        try %%% just so that whole thing does not shut down in case of temperature problem
            if param.dotemp == 1 && exist(param.tempfile) == 2
                load(seriesinfo_outputfilename)
                if ~contains(param.unique_id,'Terra')
                    %%% this is for old format of temperature files before
                    %%% Terra and when putin still showed the wrong time
                    scanner = 'Magnetom'; %make_temperature_plot(param.tempfile, seriesinfo, temperature_traceoutputfile, temperature_figureoutputfile, 'Magnetom');
                elseif contains(param.unique_id,'Terra') && str2double(param.scandate) < 200900
                    scanner = 'Terra_before_Sep20';
                else
                    %%% this is updated file
                    scanner = 'Terra'; %
                end
                try %%% there are different ways this file was formated, which changed over time, so I had to write some work arounds to accommodate this.
                    %%% will not work if the file format changes again
                    make_temperature_plot(param.tempfile, seriesinfo, temperature_traceoutputfile, temperature_figureoutputfile, scanner);
                catch
                    make_temperature_plot(param.tempfile, seriesinfo, temperature_traceoutputfile, temperature_figureoutputfile, 'Weird');
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MP2RAGE
    %%% this part i started implementing in order to use the MP2RAGE to
    %%% get an estimate of T1, which is then used in the SESTE B1+ map
    %%% calculation.
    for run = 1:length(folders.mp2rage_INV2)
        outdir = [param.processed_data_folder,'/MP2RAGE'];
        if exist(outdir) == 0
            system(['mkdir ', outdir]);
        end
        outputthere = dir([outdir,'/*',num2str(run),'*.txt']);
        if length(outputthere) == 0
            inv2_file = get_first_nifti_file_from_folder([folders.mp2rage_INV2(run).folder,'/',folders.mp2rage_INV2(run).name]);
            uni_file = get_first_nifti_file_from_folder([folders.mp2rage_UNI(run).folder,'/',folders.mp2rage_UNI(run).name]);
            if length(folders.mp2rage_T1) > 0 
                t1_file = get_first_nifti_file_from_folder([folders.mp2rage_T1(run).folder,'/',folders.mp2rage_T1(run).name]);
            else
                t1_file = '/nonsense/nonexistant.nii';
            end
            if length(dir([outdir,'/*.nii'])) == 0
                process_mp2rage(param, inv2_file, uni_file, t1_file, outdir, param.bm_thr_mp2rage, run);
            end
        end
    end     
     %%% extract value to use from last run
     if length(folders.mp2rage_INV2) ~= 0 & ~strcmp(param.useforT1ref, 'reference')
         if strcmp(param.useforT1ref, 'qmrlab')
            t1_file_to_use = [outdir,'/MP2RAGE_repeat_',sprintf('%.02d',length(folders.mp2rage_INV2)),'_T1_calculated_from_scanner_signal_masked.nii'] ;
         elseif strcmp(param.useforT1ref, 'scanner')
            t1_file_to_use = [outdir,'/MP2RAGE_repeat_',sprintf('%.02d',length(folders.mp2rage_INV2)),'_T1_calculated_with_qmrlab_signal_masked.nii'] ;
         end
         [error median] = system(['fslstats ',t1_file_to_use,' -P 50']);
          param.reference_T1 = str2double(median);
     end
             
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% B1 mapping
    %%% puts all intermediate calculation files in one folder
    b1_outputdir = [param.processed_data_folder,'/B1mapping'];
    system(['mkdir ',param.processed_data_folder]);
    system(['mkdir ',b1_outputdir]);

    %%%%%%%%%%%%% SESTE/AL B1+ %%%%%%%%%%%%%%%%%%%%%%%%%
    if find(strcmp(param.b1_method,'AL')) > 0 %%% if AL (Antoine Lutti) is specified as one of the methods in the param file, compute maps
       if length(folders.al_b1) > 0
            %%% considers all runs of the b1 acquisition
            for al_repeat = 1:length(folders.al_b1)
                %% B1 map via hMRI toolbox 
                %%% requires also B0 map
                al_outputdir=[b1_outputdir,'/B1_AL_calc_repeat_',sprintf('%02d',al_repeat)];
                curr_output_filename = [b1_outputdir,'/B1map_AL_repeat_',sprintf('%02d',al_repeat),'.nii'];
                if exist(curr_output_filename) == 0
                    system(['mkdir ',al_outputdir]);
                    curr_folders.al_b1 = folders.al_b1(al_repeat);
                    curr_folders.al_b0 = folders.al_b0(2*(al_repeat-1)+1:2*(al_repeat-1)+2); 
                    b1_files_to_use = convert_dir_output_to_cell_structure(dir([curr_folders.al_b1(end).folder,'/', curr_folders.al_b1(end).name,'/*.nii'])); %%% files from b1 mapping sequence
                    b0_magnitude_files_to_use = convert_dir_output_to_cell_structure(dir([curr_folders.al_b0(1).folder,'/', curr_folders.al_b0(1).name,'/*00001-*01-*.nii'])); %%% it created a weird other image that we ignore
                    b0_phase_files_to_use = convert_dir_output_to_cell_structure(dir([curr_folders.al_b0(2).folder,'/', curr_folders.al_b0(2).name,'/*.nii'])); %%% phase images are in second folder 
                    b0_files_to_use = [b0_magnitude_files_to_use;b0_phase_files_to_use]; %%% concatenate filenames  
                    while exist(curr_output_filename) == 0 
                        try
                           system(['yes | rm -R ',al_outputdir]); %%% just in case the folder exists, avoid having run02 etc. in there 
                        end
                        try
                            parameter_file_new = [b1_outputdir,'/AL_parameters.m'];
                            if exist(parameter_file_new) ~= 0
                                delete(parameter_file_new)
                            end
                            %%% replace assumed T1 in hmri toolbox parameter file
                            system(['sed "s/1633/',num2str(param.reference_T1),'/g" ', param.hmri_processing_parameters,' > ',parameter_file_new])
                            local_default_file = [param.bashscriptfolder,'/hmri_local_defaults.m'];
                            %%% runs b1 calculation thruogh toolbox
                            create_b1_map_with_hmri_toolbox(local_default_file,b1_files_to_use,b0_files_to_use,parameter_file_new,al_outputdir); %%% run function
                            low_res_b1_map = dir([al_outputdir,'/B1mapCalc/uB1map*.nii']); %%% b0-undistorted map without smoothing or padding
                            system(['cp ',[low_res_b1_map(1).folder,'/',low_res_b1_map(1).name],' ',curr_output_filename]);
                        end
                    end
                end
            end
       end
    end
    %%%%%%%%%%%%% AFI B1+ %%%%%%%%%%%%%%%%%%%%%%%%%
    if find(strcmp(param.b1_method,'AFI')) > 0 %%% if AFI is specified as one of the methods in the param file, compute maps
        if length(folders.afi) > 0
            for afi_repeat = 1:length(folders.afi)
                pdw_files_for_reference = convert_dir_output_to_cell_structure(dir([folders.pdw(end).folder,'/', folders.pdw(end).name,'/*.nii'])); %%% for registration purpose take first run of pdw
                afi_outputdir=[b1_outputdir,'/B1_AFI_calc_repeat_',sprintf('%02d',afi_repeat)];
                curr_output_filename = [b1_outputdir,'/B1map_AFI_repeat_',sprintf('%02d',afi_repeat),'.nii'];
                while exist(curr_output_filename) == 0
                    B1_files=convert_dir_output_to_cell_structure(dir([folders.afi(afi_repeat).folder,'/', folders.afi(afi_repeat).name,'/*.nii']));
                    try
                       system(['yes | rm -R ',afi_outputdir]); %%% just in case the folder exists, avoid having run02 etc. in there 
                    end
                    try %%% calculates but spits out error
                        display('calculate map')
                        local_default_file = [param.bashscriptfolder,'/hmri_local_defaults.m'];
                        create_afi_b1_map_with_hmri_toolbox(local_default_file, param,B1_files,pdw_files_for_reference,afi_outputdir)
                    end
                    low_res_b1_map = dir([afi_outputdir,'/B1mapCalc/*unsmoothed_B1map.nii']);
                    system(['cp ',[low_res_b1_map(1).folder,'/',low_res_b1_map(1).name],' ',curr_output_filename]);
                end
            end
        end
    end
    %%%%%%%%%%%%% MAFI B1+ %%%%%%%%%%%%%%%%%%%%%%%%%
    if find(strcmp(param.b1_method,'MAFI')) > 0 %%% this was an attempt to calculate MAFI maps which completely failed for the MALDI data
        if length(folders.mafi) > 0
            for mafi_repeat = 1:length(folders.mafi)/2
                curr_output_filename = [b1_outputdir,'/B1map_MAFI_repeat_',sprintf('%02d',mafi_repeat),'.nii'];
                while exist(curr_output_filename) == 0
                    file1 = dir([folders.mafi(mafi_repeat).folder,'/', folders.mafi(mafi_repeat).name,'/*-1.nii']); %%% may be several in folder due to some weird toolbox conversion duplication
                    file1 = [file1(1).folder,'/',file1(1).name];
                    %%% this would be as described from manual
                    %file2 = dir([folders.mafi(mafi_repeat).folder,'/', folders.mafi(mafi_repeat).name,'/*-2.nii']);
                    %file2 = [file2(1).folder,'/',file2(1).name];
                    %%% this would be taking first echo from other folder
                    file2 = dir([folders.mafi(mafi_repeat).folder,'/', folders.mafi(mafi_repeat+1).name,'/*-1.nii']);
                    file2 = [file2(1).folder,'/',file2(1).name];
                    jsonfile1 =  dir([folders.mafi(mafi_repeat).folder,'/', folders.mafi(mafi_repeat).name,'/*-1.json']);
                    [desired_fa, tr_diff] = get_mafi_parameters_from_json_file([jsonfile1(1).folder,'/',jsonfile1(1).name]);
                    create_mafi_b1_map(file2, file1, desired_fa, tr_diff, curr_output_filename)
                end
            end
        end
    end
    
    %%%%%%%%% BS B1+ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if find(strcmp(param.b1_method,'BS')) > 0 %%% this is based on something luke has done
        if length(folders.b1_bs) > 0
            for s1_bs_repeat = 1:length(folders.b1_bs)/2 %%% because always two folders
                curr_output_filename = [b1_outputdir,'/B1map_BS_repeat_',sprintf('%02d',s1_bs_repeat),'.nii'];
                curr_mask_output_filename = [b1_outputdir,'/B1map_BS_brain_mask_repeat_',sprintf('%02d',s1_bs_repeat),'.nii'];
                curr_output_filename_diff = [b1_outputdir,'/B1map_BSphasediff_repeat_',sprintf('%02d',s1_bs_repeat),'.nii'];
                first_folder_idx = 1 + (s1_bs_repeat-1) * 2;
                list_of_files_1 = dir([folders.b1_bs(first_folder_idx).folder,'/',folders.b1_bs(first_folder_idx).name,'/*.nii']);
                list_of_files_2 = dir([folders.b1_bs(first_folder_idx + 1).folder,'/',folders.b1_bs(first_folder_idx + 1).name,'/*.nii']); 
                list_of_files = [list_of_files_1; list_of_files_2];
                if exist(curr_output_filename) == 0
                    create_BS_fieldmap(list_of_files, curr_output_filename, curr_output_filename_diff, curr_mask_output_filename, 0);  %%% turn off smoothing with 0, previously 2
                end
            end
        end
    end
    
    %% copy custom b1map if a precalculated b1map was provided in the param structure
    if find(strcmp(param.b1_method,'custom')) > 0 & exist(param.custom_b1_filename) == 2
        curr_output_filename = [b1_outputdir,'/B1map_custom_repeat_01.nii'];
        system(['cp ', param.custom_b1_filename, ' ',curr_output_filename]);
    end
    
    %% make a fake B1 map if none is available
    %%% this is done so that there are no errors later, it is just a
    %%% uniform b1 of 100%
    if find(strcmp(param.b1_method,'none')) > 0
        t1_filenames = dir([folders.t1w(1).folder,'/',folders.t1w(1).name,'/*.nii']); 
        first_echo_file = [t1_filenames(1).folder,'/',t1_filenames(1).name];
        curr_output_filename = [b1_outputdir,'/B1map_none_repeat_01.nii'];
        if exist(curr_output_filename) == 0
            system(['fslmaths ', first_echo_file, ' -abs -bin -mul 100 ', curr_output_filename]);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% calculate MPMs for each run separately
    %%% these are all stored in the MPM_calc directory with file names that
    %%% do indicate which processing and corrections have been done
    mpm_output_dir = [param.processed_data_folder,'/MPM_calc'];
    system(['mkdir ', mpm_output_dir]);
    runs_to_consider = 1:length(folders.t1w);

    %% start MPM processing
    for mpm_run = runs_to_consider

            %% create brain mask and capsule mask based on a high intensity file
            signal_mask_file = [mpm_output_dir,'/signal_mask_',param.resolution,'_run',sprintf('%02d',mpm_run),'.nii']; %%% this is brain mask + water capsule (if capsule was measured)
            brain_mask_file = [mpm_output_dir,'/brain_mask_',param.resolution,'_run',sprintf('%02d',mpm_run),'.nii'];
            if exist(param.own_brain_mask_file) == 2 %%% if it is specified in param that an already existing brain mask file should be used
                system(['cp ',param.own_brain_mask_file,' ',brain_mask_file]);
                system(['cp ',param.own_brain_mask_file,' ',signal_mask_file]);
            else %%% make brain mask using the specified threshold
                brain_mask_file_diluted = [mpm_output_dir,'/brain_mask_',param.resolution,'_run',sprintf('%02d',mpm_run),'_diluted.nii'];
                capsule_mask_file = [mpm_output_dir,'/capsule_mask_',param.resolution,'_run',sprintf('%02d',mpm_run),'.nii']; %%% this may end up empty
                t1_filenames = dir([folders.t1w(mpm_run).folder,'/',folders.t1w(mpm_run).name,'/*.nii']); 
                first_echo_file = [t1_filenames(1).folder,'/',t1_filenames(1).name];
                if exist(brain_mask_file) == 0 || exist(signal_mask_file) == 0
                    %%% threshold the first echo of t1w with the signal
                    %%% threshold specified in param
                    system(['fslmaths ',first_echo_file,' -thr ', num2str(param.bm_thr), ' -bin ',signal_mask_file]); %%% threshold based on signal intensity of raw file
                    system(['bash ', param.bashscriptfolder, '/take_largest_cluster_from_a_mask ',signal_mask_file,' ',brain_mask_file, ' ', mpm_output_dir]); %%% get rid of noisy voxels around          
                    if param.mtv_calc == 1
                        %%% make capsule mask by trying to grab second largest
                        %%% cluster
                        system(['fslmaths ', signal_mask_file,' -sub ', brain_mask_file,' ',capsule_mask_file]);
                        system(['bash ', param.bashscriptfolder, '/take_largest_cluster_from_a_mask ',capsule_mask_file,' ',capsule_mask_file, ' ', mpm_output_dir]); %%% get rid of noisy voxels around
                        system(['fslmaths ', brain_mask_file, ' -add ', capsule_mask_file, ' -bin ', signal_mask_file]);
                    else
                        system(['cp ', brain_mask_file, ' ', signal_mask_file]);
                    end
                    %%% make diluted brain mask to fill some holes
                    if exist(brain_mask_file_diluted) == 0
                        system(['bash ',param.bashscriptfolder,'/fill_gaps_in_brain_mask ',brain_mask_file,' ',brain_mask_file_diluted,' 5 ', mpm_output_dir]);
                    end
                end
            end
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% define how the file names should be called, depending on
            %%% what processing was specified in param (this is useful in
            %%% case a loop is run over parameters to compare the effect of
            %%% different steps)
            processed_file_string = {'processed_s'};
%             if param.do_denoising == 1
%                 processed_file_string = strcat(processed_file_string,'_dn');
%             end
            if param.do_movement_correction == 1 %%% in case sample moved between sequences
                processed_file_string = strcat(processed_file_string,'_reg');
            end
            if param.do_receive_bias_correction == 1 %%% just relevant for PD
                processed_file_string = strcat(processed_file_string,'_rbc');
            end
            if param.do_distortion_correction == 1 %%% to account for tiny shifts between odd and even echos
                processed_file_string = strcat(processed_file_string,'_u');
            end
            processed_file_string =  processed_file_string{1};
            processed_file_folder = [param.processed_data_folder,'/', processed_file_string];
            if exist(processed_file_folder) == 0
                system(['mkdir ', processed_file_folder]);
            end
            %% copy all raw files to the processed directory
            try %%% if we have mtw
                 all_folders = {folders.t1w(mpm_run).name, folders.pdw(mpm_run).name, folders.mtw(mpm_run).name};
            catch
                 all_folders = {folders.t1w(mpm_run).name, folders.pdw(mpm_run).name};
            end
            if length(folders.additionalFA) > 0 %%% if another weighting was acquired (usually that would be ernst angle)
                try
                    all_folders = {folders.pdw(mpm_run).name, folders.t1w(mpm_run).name, folders.mtw(mpm_run).name, folders.additionalFA(mpm_run).name};
                catch
                    all_folders = {folders.pdw(mpm_run).name, folders.t1w(mpm_run).name, folders.additionalFA(mpm_run).name};
                end
            end
            for folder = 1:length(all_folders)
                if exist([processed_file_folder,'/',all_folders{folder}]) == 0
                    system(['cp -R ', param.raw_data_folder_nifti,'/',all_folders{folder}, ' ', processed_file_folder]);
                end
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% movement correction
            if param.do_movement_correction == 1
                    for folder = 1:length(all_folders)
                      to_corr_foldername = [processed_file_folder,'/',all_folders{folder}];
                      system(['mkdir ', to_corr_foldername])
                      to_corr_filenames = dir([to_corr_foldername,'/s*.nii']);
                      to_corr_jsonfilenames = dir([to_corr_foldername,'/s*.json']);
                      if folder == 1 %%% this is the first sequence which is used as reference sequence
                          %%% in the param, it is specified which echo
                          %%% should be the reference echo
                          reference_image = [to_corr_filenames(param.echo2reg).folder,'/',to_corr_filenames(param.echo2reg).name]; %%% change referene image to be a file with more
                          %%% contrast!
                          %%% just copy files
                          for f = 1:length(to_corr_filenames)
                            niftifilename_in = to_corr_filenames(f).name;
                            niftifilename_out = ['reg_', niftifilename_in];
                            if exist([to_corr_foldername,'/',niftifilename_out]) == 0 && strcmp(niftifilename_in(1),'s') == 1
                                system(['cp ', to_corr_foldername,'/',niftifilename_in, ' ', to_corr_foldername,'/',niftifilename_out]);
                            end
                            jsonfilename_in = to_corr_jsonfilenames(f).name;
                            jsonfilename_out = ['reg_', jsonfilename_in];
                            if exist(jsonfilename_out) == 0 && strcmp(niftifilename_in(1),'s') == 1
                                system(['cp ', to_corr_foldername,'/',jsonfilename_in, ' ', to_corr_foldername,'/',jsonfilename_out]);
                            end
                          end
                      else %%% now the images in the other folders are registered
                          registration_img = [to_corr_filenames(param.echo2reg).folder,'/',to_corr_filenames(param.echo2reg).name];
                          reg_matrix = [to_corr_foldername,'/affine'];
                          %%% use flirt to compute a registration matrix
                          %%% with the parameter file specified DOF
                          if exist(reg_matrix) == 0
                            system(['flirt -in ', registration_img, ' -ref ', reference_image ,' -dof ',num2str(param.mocodof),' -omat ', reg_matrix]);
                          end
                          for f = 1:length(to_corr_filenames)
                            niftifilename_in = to_corr_filenames(f).name;
                            niftifilename_out = ['reg_', niftifilename_in];
                            %%% apply the registration matrix
                            if exist([to_corr_foldername,'/',niftifilename_out]) == 0 && strcmp(niftifilename_in(1),'s') == 1
                                system(['flirt -in ', to_corr_foldername,'/',niftifilename_in, ' -ref ', reference_image, ' -applyxfm -init ', reg_matrix , ' -out ', to_corr_foldername,'/',niftifilename_out]);
                            end
                            jsonfilename_in = to_corr_jsonfilenames(f).name;
                            jsonfilename_out = ['reg_', jsonfilename_in];
                            %%% json file just copied over
                            if exist(jsonfilename_out) == 0 && strcmp(niftifilename_in(1),'s') == 1
                                system(['cp ', to_corr_foldername,'/',jsonfilename_in, ' ', to_corr_foldername,'/',jsonfilename_out]);
                            end
                          end
                      end
                    end

             end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% receive bias correction
            %%% 1: make highres rb map
            %%% this is based on the comparison of low resolution images
            %%% acquired with the receive vs transmit coil
            %%% for the Terra scanner some of the echos in these sequences
            %%% were sometimes very corrupted so this does not work as well
            %%% anymore. 
            if param.do_receive_bias_correction == 1 && length(folders.receive_bias_receive) > 0 && length(folders.receive_bias_transmit) > 0 
                receive_bias_map_to_use = [mpm_output_dir,'/receive_bias_average_map.nii'];
                high_res_rb_map = [mpm_output_dir,'/receive_bias_map_highres_',param.resolution,'_run',sprintf('%02d',mpm_run),'.nii']; %%% decide on output file name
                if exist(receive_bias_map_to_use) == 0 
                    for repeat = 1:length(folders.receive_bias_receive) %%% consider various runs
                        receive_nifti_folder = [folders.receive_bias_receive(repeat).folder,'/',folders.receive_bias_receive(repeat).name];
                        transmit_nifti_folder = [folders.receive_bias_transmit(repeat).folder,'/',folders.receive_bias_transmit(repeat).name];
                        output_filename = [mpm_output_dir,'/receive_bias_repeat_',sprintf('%02d',repeat),'.nii'];
                        if exist(output_filename) == 0
                            output_filename_for_trouble_shoot = [mpm_output_dir,'/receive_bias_repeat_',sprintf('%02d',repeat),'_allechos.nii'];
                            temp_folder = mpm_output_dir;
                            calculate_receive_bias_map(transmit_nifti_folder, receive_nifti_folder, output_filename);
                          %  system(['bash ', param.bashscriptfolder, '/calculate_receive_bias_map_without_smooth ',receive_nifti_folder,' ',transmit_nifti_folder,' ',temp_folder,' ',output_filename,' ', output_filename_for_trouble_shoot]);
                            end
                    end
                    %%% average runs
                    system(['fslmerge -t ',mpm_output_dir,'/receive_bias_all_repeats.nii ',mpm_output_dir,'/receive_bias_rep*.nii']);
                    system(['fslmaths ',mpm_output_dir,'/receive_bias_all_repeats.nii -Tmean ',receive_bias_map_to_use]);
                end
                if exist(high_res_rb_map) == 0 
                    %%% process this map like the b1map
                    system(['bash ', param.bashscriptfolder, '/new_b1_processing_median_filter ',receive_bias_map_to_use,' ',signal_mask_file,' ',high_res_rb_map,' ',param.rbmap_smooth_in_mm_median,' ',param.rbmap_smooth_in_mm,' ',mpm_output_dir]);
                end
                
                %%% 2: correct data (for scripting simplicity the weighted files are
                %%% corrected before MPMs calculated)
                for folder = 1:length(all_folders)
                      to_corr_foldername = [processed_file_folder,'/',all_folders{folder}];
                      if length(dir([to_corr_foldername, '/u2*2.*'])) > 0
                        system(['yes | rm ', to_corr_foldername, '/u2*2.*'])
                      end
                      %%% define what the input filenames should be
                      if param.do_movement_correction == 1
                         to_corr_filenames = dir([to_corr_foldername,'/reg*.nii']);
                         to_corr_jsonfilenames = dir([to_corr_foldername,'/reg*.json']);
                      else
                         to_corr_filenames = dir([to_corr_foldername,'/s*.nii']);
                         to_corr_jsonfilenames = dir([to_corr_foldername,'/s*.json']);
                      end
                      for f = 1:length(to_corr_filenames)
                        niftifilename_in = to_corr_filenames(f).name;
                        niftifilename_out = ['rbc_', niftifilename_in];
                        if exist([to_corr_foldername,'/',niftifilename_out]) == 0 && (strcmp(niftifilename_in(1),'s') == 1 || strcmp(niftifilename_in(1),'r') == 1)
                            %%% correct by dividing weighted images by
                            %%% receive bias map
                            system(['fslmaths ', to_corr_foldername,'/',niftifilename_in, ' -div ', high_res_rb_map, ' ', to_corr_foldername,'/',niftifilename_out]);
                        end
                        jsonfilename_in = to_corr_jsonfilenames(f).name;
                        jsonfilename_out = ['rbc_', jsonfilename_in];
                        if exist(jsonfilename_out) == 0 && (strcmp(niftifilename_in(1),'s') == 1 || strcmp(niftifilename_in(1),'r') == 1)
                            %%% again just copy over the json file so that
                            %%% it exists
                            system(['cp ', to_corr_foldername,'/',jsonfilename_in, ' ', to_corr_foldername,'/',jsonfilename_out]);
                        end
                      end
                end
            end
            

            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% distortion correction
            if param.do_distortion_correction == 1
                for folder = 1:length(all_folders)
                      to_corr_foldername = [processed_file_folder,'/',all_folders{folder}];
                      %%% complicated way to define what the input filename
                      %%% should be called like
                      if param.do_receive_bias_correction == 1
                             to_corr_filenames = dir([to_corr_foldername,'/rbc*.nii']);
                             to_corr_jsonfilenames = dir([to_corr_foldername,'/rbc*.json']);
                      else
                          if param.do_movement_correction == 1
                            to_corr_filenames = dir([to_corr_foldername,'/reg*.nii']);
                                to_corr_jsonfilenames = dir([to_corr_foldername,'/reg*.json']);
                          else
                             to_corr_filenames = dir([to_corr_foldername,'/s*.nii']);
                             to_corr_jsonfilenames = dir([to_corr_foldername,'/s*.json']);
                          end
                      end
                      output_geom = [processed_file_folder,'/Geom_',all_folders{folder},'.nii'];
                      output_bfield = [processed_file_folder,'/HySCOv2_Geom_',all_folders{folder},'.nii'];
                      check_if_output_exists = dir([to_corr_foldername,'/u*.nii']);
                      if length(check_if_output_exists) < length(to_corr_filenames)
                          %%% use nik's script to perform the correction
                        apply_niks_distortion_correction(convert_dir_output_to_cell_structure(to_corr_filenames), output_geom, output_bfield, to_corr_foldername, param.permdim)
                      end
                      %%% copy over json filenames so that corrected files
                      %%% can be used
                      for f = 1:length(to_corr_jsonfilenames)
                         jsonfilename_in = to_corr_jsonfilenames(f).name;
                         jsonfilename_out = ['u', jsonfilename_in];
                         if exist([to_corr_foldername,'/',jsonfilename_out]) == 0 && (strcmp(jsonfilename_in(1),'u') == 0)
                            system(['cp ', to_corr_foldername,'/',jsonfilename_in, ' ', to_corr_foldername,'/',jsonfilename_out]);
                         end
                      end
                end
            end

            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% now that we have a brainmask, process the B1 maps properly, we are already in correct param.resolution
            for b = 1:length(param.b1_method)
                b1_map_type = param.b1_method{b};
                all_raw_maps = dir([b1_outputdir,'/B1map_',b1_map_type,'_repeat*.nii']);
                if length(all_raw_maps) > 0 && exist([b1_outputdir,'/B1map_',b1_map_type,'_processed_',param.resolution,'_mean.nii']) == 0
                    %%% for all repeats for now to test
                    for repeat = 1:length(all_raw_maps)
                        map_to_process = [all_raw_maps(repeat).folder,'/',all_raw_maps(repeat).name];
                        new_map = [b1_outputdir,'/B1map_',b1_map_type,'_processed_',param.resolution,'_repeat',sprintf('%02d',repeat),'.nii']
                        %padding_parameter_in_vxs = 3; %%% number of voxels being padded roughly
                        %process_b1_map(map_to_process, signal_mask_file, smooth_in_mm, padding_parameter_in_vxs, new_map);
                        if strcmp(param.b1_method_to_use,'none')
                            %%% the uniform b1map only has to be upsampled
                            system(['flirt -ref ', signal_mask_file, ' -in ', map_to_process, ' -out ', new_map , ' -applyxfm -usesqform -interp nearestneighbour ']);                        
                        else %%% for all actual b1 maps an optimised processing is applied (to smooth across the very strong boundaries)
                           % ['bash ', param.bashscriptfolder, '/new_b1_processing ',map_to_process,' ',signal_mask_file,' ',new_map,' ',num2str(param.b1map_smooth_in_mm),' ',b1_outputdir]
                            system(['bash ', param.bashscriptfolder, '/new_b1_processing ',map_to_process,' ',signal_mask_file,' ',new_map,' ',num2str(param.b1map_smooth_in_mm),' ',b1_outputdir]);
                        end
                    end
                    %%% create average map and CVs (coefficient of variation maps across the repetition) for unprocessed
                    system(['fslmerge -t ',b1_outputdir,'/B1map_',b1_map_type,'_all_repeats ',b1_outputdir,'/B1map_',b1_map_type,'_repeat*.nii']);
                    system(['fslmaths ',b1_outputdir,'/B1map_',b1_map_type,'_all_repeats.nii -Tmean ',b1_outputdir,'/B1map_',b1_map_type,'_mean.nii']);
                    system(['fslmaths ',b1_outputdir,'/B1map_',b1_map_type,'_all_repeats.nii -Tstd ',b1_outputdir,'/B1map_',b1_map_type,'_std.nii']);
                    system(['fslmaths ',b1_outputdir,'/B1map_',b1_map_type,'_std.nii -div ',b1_outputdir,'/B1map_',b1_map_type,'_mean.nii -mul 100 ',b1_outputdir,'/B1map_',b1_map_type,'_CVmap.nii']);
                    %process_b1_map([b1_outputdir,'/B1map_',b1_map_type,'_mean.nii'], signal_mask_file, smooth_in_mm, [b1_outputdir,'/B1map_',b1_map_type,'_processed_',param.resolution,'_mean.nii']);
                    %%% create average map and CVs for processed
                    system(['fslmerge -t ',b1_outputdir,'/B1map_',b1_map_type,'_processed_',param.resolution,'_all_repeats ',b1_outputdir,'/B1map_',b1_map_type,'_processed_',param.resolution,'_repeat*.nii']);
                    system(['fslmaths ',b1_outputdir,'/B1map_',b1_map_type,'_processed_',param.resolution,'_all_repeats.nii -Tmean ',b1_outputdir,'/B1map_',b1_map_type,'_processed_',param.resolution,'_mean.nii']);
                    system(['fslmaths ',b1_outputdir,'/B1map_',b1_map_type,'_processed_',param.resolution,'_all_repeats.nii -Tstd ',b1_outputdir,'/B1map_',b1_map_type,'_processed_',param.resolution,'_std.nii']);
                    system(['fslmaths ',b1_outputdir,'/B1map_',b1_map_type,'_processed_',param.resolution,'_std.nii -div ',b1_outputdir,'/B1map_',b1_map_type,'_processed_',param.resolution,'_mean.nii -mul 100 ',b1_outputdir,'/B1map_',b1_map_type,'_processed_',param.resolution,'_CVmap.nii']);
                    
                end
            end

            %% IMPLEMENT B1 COMPARISON HERE (in case several b1mapping approaches were used, a figure is created to illustrate correspondence)
            b1_comp_figure_filename = [b1_outputdir,'/',param.unique_id,'_B1map_comparison_',param.resolution,'.png'];
            if length(param.b1_method)>1 && exist(b1_comp_figure_filename) == 0
                subpl_count = 1;
                for b1 = 1:length(param.b1_method)
                    b1_map1_filename = [b1_outputdir,'/B1map_',param.b1_method{b1},'_processed_',param.resolution,'_mean.nii'];
                    for b2 = 1:length(param.b1_method)
                        b1_map2_filename = [b1_outputdir,'/B1map_',param.b1_method{b2},'_processed_',param.resolution,'_mean.nii'];
                        %%% shuff into function
                        f = figure(99);
                        set(f,'position',[0 0 2000 800]);
                        set(gcf,'Visible','off');
                        subplot(length(param.b1_method),length(param.b1_method),subpl_count);
                        subpl_count = subpl_count + 1;
                        compare_B1_maps_with_BA(b1_map1_filename, b1_map2_filename);
                        xlabel(['mean (',param.b1_method{b1},'&',param.b1_method{b2},')'])
                        ylabel(['diff (',param.b1_method{b1},'-',param.b1_method{b2},')'])       
                    end
                end
                try
                    saveas(f,b1_comp_figure_filename); 
                end
            end
           
            
            %% COMPUTE MPMS
          
            if length(folders.t1w) > 0 && length(folders.pdw) > 0 %%% & length(param.mt_folder_series_numbers) > 0 
        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% read B1 map
                high_res_b1_map = [b1_outputdir,'/B1map_',param.b1_method_to_use,'_processed_',param.resolution,'_mean.nii'];
                B1vol = spm_vol(high_res_b1_map);
                percB1 = spm_read_vols(B1vol);
               
                %% now calculate MPMs 
                %%% but only if they do not exist yet
                if 1 == 1 %exist([mpm_output_dir,'/MTsat_real_TE0_WLS_',processed_file_string,'_',param.resolution,'_run',sprintf('%02d',mpm_run),'_brain_masked.nii']) == 0 %% && exist([mpm_output_dir,'/R1_rbc_',param.resolution,'_run',sprintf('%02d',mpm_run),'_brain_masked.nii']) == 0

                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% create data frame in the right format
                    %%% from processed file string take last string between _
                    %%% and use this as prefix for filename

                    proc_parts = strsplit(processed_file_string, '_'); %%% all elements of processing
                
                    for echo_iteration = 1:length(param.echos_to_use) %%% this specifies which echos should be used for calculation. skipfirst was introduced for saturated images
                       whatecho = param.echos_to_use{echo_iteration}; %%% either all or number of echos to use
                       if ~strcmp(whatecho, 'all') & ~strcmp(whatecho, 'skipfirst')
                            processed_file_string_with_echos = [processed_file_string,'_',sprintf('%.02d',whatecho),'echos'];
                       else
                           processed_file_string_with_echos = processed_file_string; %%% all does not need to be specified
                       end
                       
                        %% find filenames and make a useable data structure for each weighting
                        weighted_fieldnames = {'t1w', 'pdw','mtw','additionalFA'}
                        for f = 1:length(weighted_fieldnames);
                            weighted_fieldname = weighted_fieldnames{f};
                            try
                                weighted.filenames.(weighted_fieldname) = dir([processed_file_folder,'/',folders.(weighted_fieldname)(mpm_run).name,'/',proc_parts{end},'*.nii']);
                                weighted.jsonfilenames.(weighted_fieldname) = dir([processed_file_folder,'/',folders.(weighted_fieldname)(mpm_run).name,'/',proc_parts{end},'*.json']);
                                if strcmp(whatecho, 'all')
                                    weighted.data.(weighted_fieldname) = make_mpm_file_structure_from_filenames(weighted.filenames.(weighted_fieldname));
                                    first_echo_file = [weighted.jsonfilenames.(weighted_fieldname)(1).folder, '/', weighted.jsonfilenames.(weighted_fieldname)(1).name];
                                elseif strcmp(whatecho, 'skipfirst')
                                    weighted.data.(weighted_fieldname) = make_mpm_file_structure_from_filenames(weighted.filenames.(weighted_fieldname)(2:end));
                                    weighted.filenames.(weighted_fieldname) = weighted.filenames.(weighted_fieldname)(2:end);
                                    weighted.jsonfilenames.(weighted_fieldname) =  weighted.jsonfilenames.(weighted_fieldname)(2:end);
                                    first_echo_file = [weighted.jsonfilenames.(weighted_fieldname)(2).folder, '/', weighted.jsonfilenames.(weighted_fieldname)(2).name];
                                else %% 1:number specified
                                    weighted.data.(weighted_fieldname) = make_mpm_file_structure_from_filenames(weighted.filenames.(weighted_fieldname)(1:whatecho));
                                    weighted.filenames.(weighted_fieldname) = weighted.filenames.(weighted_fieldname)(1:whatecho);
                                    weighted.jsonfilenames.(weighted_fieldname) =  weighted.jsonfilenames.(weighted_fieldname)(1:whatecho);
                                    first_echo_file = [weighted.jsonfilenames.(weighted_fieldname)(1).folder, '/', weighted.jsonfilenames.(weighted_fieldname)(1).name];
                                end
                                
                                [weighted.data.(weighted_fieldname).desired_fa, weighted.data.(weighted_fieldname).sinc_pulse_duration, weighted.data.(weighted_fieldname).srf01, weighted.data.(weighted_fieldname).srf02] = get_fa_parameters_from_json_file(first_echo_file);
                            end
                        end
                        
                        if param.correct_fa == 1 %%% this was introduced to check the voltages in order to correct the flip angle (in case the desired flip angle - that is noted in the json file - was not actually achieved)
                            effective_t1w_fa = weighted.data.pdw.desired_fa * weighted.data.t1w.srf01/weighted.data.pdw.srf01 * weighted.data.t1w.sinc_pulse_duration/weighted.data.pdw.sinc_pulse_duration;
                        else
                            effective_t1w_fa = weighted.data.t1w.desired_fa;
                        end
                        if ~isnan(effective_t1w_fa) %%% write the actual flip angle out
                            dlmwrite([mpm_output_dir,'/effective_t1w_fa_',param.resolution,'_run',sprintf('%02d',mpm_run)],effective_t1w_fa)
                        else
                            effective_t1w_fa = weighted.data.t1w.desired_fa;
                        end
                        weighted.data.t1w.fa = deg2rad(effective_t1w_fa);

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %% calculate R2s with toolbox

                        TR = weighted.data.t1w(1).TR; % assumed equal for all acquisition
                        try
                            MAPS.T1w_echo8 = weighted.data.t1w.data(:,:,:,8);
                        end
                        %%% first implement an OLS fit. if mask for bubbles
                        %%% is desired then the model fit will also be
                        %%% provided. this is computationally expensive
                        %%% (due to inefficient implementation), so only
                        %%% done when explicitely desired
                        if param.createbubblemask == 1
                            try
                                [R2s_OLS,varexpl,PDw0_OLS,T1w0_OLS,MTw0_OLS] = weighted2R2s_IL_model_fit(weighted.data.pdw,weighted.data.t1w,weighted.data.mtw); 
                            catch
                                [R2s_OLS,varexpl,PDw0_OLS,T1w0_OLS] = weighted2R2s_IL_model_fit(weighted.data.pdw,weighted.data.t1w); 
                            end
                            R2s_OLS.data(isinf(R2s_OLS.data)) = 0;
                            R2s_OLSvals = R2s_OLS.data(:);
                            idx = find(R2s_OLSvals > 0);
                            thresh = nanmean(R2s_OLSvals) + 2*nanstd(R2s_OLSvals); %%% establish that values in suspected bubbles need to be outliers in distribution
                            %%% create a preliminary mask by using voxels
                            %%% with less than 98% of variance explained
                            %%% and R2* values that are more than 2 std
                            %%% higher than the mean
                            MAPS.bubble_mask_preliminary = varexpl < .98 & R2s_OLS.data > thresh;
                            MAPS.OLSfit_R2 = varexpl;
                        else
                            try
                                [R2s_OLS,PDw0_OLS,T1w0_OLS,MTw0_OLS] = weighted2R2s(weighted.data.pdw,weighted.data.t1w,weighted.data.mtw); 
                            catch
                                [R2s_OLS,PDw0_OLS,T1w0_OLS] = weighted2R2s(weighted.data.pdw,weighted.data.t1w); 
                            end
                        end
                        MAPS.R2s_OLS = R2s_OLS.data * 1000; %%% convert to ms
                        MAPS.PDw0_OLS = PDw0_OLS.data;
                        MAPS.T1w0_OLS = T1w0_OLS.data;

                        %% this is not used routinely but can be adapted if separate R2s for the different contrasts are desired
                         %%% calculate R2s for each contrast with hMRI
                        %%% toolbox
                         wims = fieldnames(weighted.filenames);
                         if 1 == 2
                         for i = 1:length(wims)
                             con = wims{i};
                             hmri_output_dir = [mpm_output_dir,'/hmri_toolbox_folder_',processed_file_string_with_echos,'_',param.resolution,'_run',sprintf('%02d',mpm_run),'_only_',con];  
                            local_default_file = [param.bashscriptfolder,'/hmri_local_defaults.m'];
                            if exist([hmri_output_dir,'/Results']) == 0
                                MPMcalc_basic_with_hmri_toolbox(local_default_file, hmri_output_dir, convert_dir_output_to_cell_structure(weighted.filenames.(con)), convert_dir_output_to_cell_structure(weighted.filenames.(con)), convert_dir_output_to_cell_structure(weighted.filenames.(con))); %%% trick the toolbox and provide twice the same
                            end
                             r2s_fname = dir([hmri_output_dir,'/Results/*R2s_WOLS.nii']);
                            try
                                mapname = ['R2s_WLS_only_',con];                          
                                MAPS.(mapname) = spm_read_vols(spm_vol([r2s_fname.folder,'/',r2s_fname.name]));
                            end
                        end
                         end
                            
                        
                        %% test luke's own R2s implementation
%                         [testr2s extrap] = hmri_calc_R2s_luke({weighted.data.pdw, weighted.data.t1w});
%                         %%% change data structure so that it fits with
%                         %%% luke's new function
%                         weighted.data.t1w2 = weighted.data.pdw; weighted.data.t1w2.TE = weighted.data.t1w2.TEs;
%                         weighted.data.pdw2 = weighted.data.pdw; weighted.data.pdw2.TE = weighted.data.pdw2.TEs;
%                         a(1) =  weighted.data.t1w2; a(2) =  weighted.data.pdw2;
%                         [testr2s extrap] = hmri_calc_R2s_luke(a,'wls3');
%                         MAPS.R2sLukeTest = testr2s;
%                         [testr2s extrap] = hmri_calc_R2s_luke(a,'ols');
%                         MAPS.R2sLukeTestOls = testr2s;
                        
                        %% calculate R2s with hMRI toolbox
                         hmri_output_dir = [mpm_output_dir,'/hmri_toolbox_folder_',processed_file_string_with_echos,'_',param.resolution,'_run',sprintf('%02d',mpm_run)];  
                        local_default_file = [param.bashscriptfolder,'/hmri_local_defaults.m'];

                        %system(['yes | rm -R ', hmri_output_dir]); %%% rerun calculations
                        %%% run pipeline to get output folder
                        if exist([hmri_output_dir,'/Results']) == 0
                            
                            if isfield(weighted.filenames,'mtw') ~= 0
                                %MPMcalc_basic_with_hmri_toolbox(local_default_file, hmri_output_dir, convert_dir_output_to_cell_structure(weighted.filenames.mtw), convert_dir_output_to_cell_structure(weighted.filenames.pdw), convert_dir_output_to_cell_structure(weighted.filenames.t1w))
                                MPMcalc_basic_with_hmri_toolbox_b1corr(local_default_file, hmri_output_dir, convert_dir_output_to_cell_structure(weighted.filenames.mtw), convert_dir_output_to_cell_structure(weighted.filenames.pdw), convert_dir_output_to_cell_structure(weighted.filenames.t1w),high_res_b1_map)
                            else
                                %MPMcalc_basic_with_hmri_toolbox(local_default_file, hmri_output_dir, {}, convert_dir_output_to_cell_structure(weighted.filenames.pdw), convert_dir_output_to_cell_structure(weighted.filenames.t1w))
                                MPMcalc_basic_with_hmri_toolbox_b1corr(local_default_file, hmri_output_dir, {}, convert_dir_output_to_cell_structure(weighted.filenames.pdw), convert_dir_output_to_cell_structure(weighted.filenames.t1w),high_res_b1_map)
                            end
                        end

                        %%% specifiy filenames to take from output folder
                         r2s_fname = dir([hmri_output_dir,'/Results/*R2s_WOLS.nii']);
                        pdw0_fname = dir([hmri_output_dir,'/Results/Supplementary/*PDw_WOLSfit_TEzero.nii']);
                        t1w0_fname = dir([hmri_output_dir,'/Results/Supplementary/*T1w_WOLSfit_TEzero.nii']);
                        mtw0_fname = dir([hmri_output_dir,'/Results/Supplementary/*MTw_WOLSfit_TEzero.nii']);

                        
                        %%% read in the files
                        
                         MPMvol = spm_vol([r2s_fname.folder,'/',r2s_fname.name]);

                        R2s_WLS = R2s_OLS; %%% copy info 
                        R2s_WLS.data = spm_read_vols(spm_vol([r2s_fname.folder,'/',r2s_fname.name]));

                        PDw0_WLS = PDw0_OLS;
                        PDw0_WLS.data = spm_read_vols(spm_vol([pdw0_fname.folder,'/',pdw0_fname.name]));

                        T1w0_WLS = T1w0_OLS;
                        T1w0_WLS.data = spm_read_vols(spm_vol([t1w0_fname.folder,'/',t1w0_fname.name]));

                        try
                            MTw0_WLS = PDw0_OLS;
                            MTw0_WLS.data = spm_read_vols(spm_vol([mtw0_fname.folder,'/',mtw0_fname.name]));
                        end

                        MAPS.R2s_WLS = R2s_WLS.data;
                        MAPS.PDw0_WLS = PDw0_WLS.data;
                        MAPS.T1w0_WLS = T1w0_WLS.data;

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %% calculate A0
                        %%% takes PDw0 and T1w0 as calculated above with the
                        %%% respetive echos
                        A0_OLS = PDwT1w2A(PDw0_OLS,T1w0_OLS,'exact',percB1/100); %%% B1 needs to be in fractions not percent
                        A0_WLS = PDwT1w2A(PDw0_WLS,T1w0_WLS,'exact',percB1/100); %%% B1 needs to be in fractions not percent
                        MAPS.AO_OLS = A0_OLS.data;
                        MAPS.A0_WLS = A0_WLS.data;

                        %% use A0 to calculate MTV
                        if param.mtv_calc == 1
                            %%% get value of pure water (with useful relaxation time) by
                            %%% getting robust estimate from within capsule
                            capsule_mask = spm_read_vols(spm_vol(capsule_mask_file));
                            capsule_mask = capsule_mask(:);
                            capsule_idx = find(capsule_mask == 1);
                            %%% OLS:
                            capsule_water_OLS = prctile(A0_OLS.data(capsule_idx),50);
                            pure_water_OLS = capsule_water_OLS / param.mtv_refdensity;
                            MAPS.MTV_OLS = 1 - A0_OLS.data / pure_water_OLS;
                            %%% WLS:
                            capsule_water_WLS = prctile(A0_WLS.data(capsule_idx),50);
                            pure_water_WLS = capsule_water_WLS / param.mtv_refdensity;
                            MAPS.MTV_OLS = 1 - A0_WLS.data / pure_water_WLS;
                        end

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %% calculate R1, either just using the first echo (echo1) or the interpolation to te = 0 (TE0, either obtained from OLS or from WLS)
                        PDw_echo1 = averageEchoes(weighted.data.pdw,1);
                        T1w_echo1 = averageEchoes(weighted.data.t1w,1);
                        try
                          additionalFAw_echo1 = averageEchoes(weighted.data.additionalFA,1);  
                        end
                        if length(folders.additionalFA) > 0 %%% if an additional weighitng was used
                            %additionalFA_echo1 = averageEchoes(weighted.data.additionalFAw,1);
                            %[Aavg,R1_not_used] = PDwT1w2quant(PDw_echo1,T1w_echo1,TR,'exact',percB1/100); %%% Luke's function requires the B1 map to be in fraction, not in percent
                            [A_echo1_additionalFA,R1_echo1_additionalFA] = weighted2AR1({PDw_echo1,T1w_echo1,additionalFAw_echo1},TR,'exact',percB1/100); 
                            %%% also produce an B1 uncorrected map as comparison
                            [A_echo1_b1uncorr_additionalFA,R1_echo1_additionalFA_b1uncorr] = weighted2AR1({PDw_echo1,T1w_echo1,additionalFAw_echo1},TR,'exact',percB1.^0); %%% shuff map of ones into b1 spot 
                        end
                        %[Aavg,R1_not_used] = PDwT1w2quant(PDw_echo1,T1w_echo1,TR,'exact',percB1/100); %%% compute T1 now; Luke's function requires the B1 map to be in fraction, not in percent
                        [A_echo1,R1_echo1] = weighted2AR1({PDw_echo1,T1w_echo1},TR,'exact',percB1/100);
                        [A_echo1_b1uncorr,R1_echo1_b1uncorr] = weighted2AR1({PDw_echo1,T1w_echo1},TR,'exact',percB1.^0);
                        [A_TE0_WLS,R1_TE0_WLS] = weighted2AR1({PDw0_WLS,T1w0_WLS},TR,'exact',percB1/100);
                        [A_TE0_WLS_b1uncorr,R1_TE0_WLS_b1uncorr] = weighted2AR1({PDw0_WLS,T1w0_WLS},TR,'exact',percB1.^0);
                        [A_TE0_OLS,R1_TE0_OLS] = weighted2AR1({PDw0_OLS,T1w0_OLS},TR,'exact',percB1/100);

                        %%% save the MPMs in a MAP structure that has all
                        %%% of them
                        MAPS.A_echo1 = A_echo1.data;
                        MAPS.R1_echo1 = R1_echo1.data;
                        MAPS.A_echo1_b1uncorr = A_echo1_b1uncorr.data;
                        MAPS.R1_echo1_b1uncorr = R1_echo1_b1uncorr.data;
                        MAPS.R1_TE0_WLS = R1_TE0_WLS.data;
                        MAPS.R1_TE0_OLS = R1_TE0_OLS.data;
                        MAPS.A_TE0_WLS_b1uncorr = A_TE0_WLS_b1uncorr.data;
                        MAPS.R1_TE0_WLS_b1uncorr = R1_TE0_WLS_b1uncorr.data;

                        %% if specified, then do slice profile correction
                        if (strcmp(param.resolution,'0p4') || strcmp(param.resolution,'0p3')) && param.dosliceprofilecorrection == 1
                            %%% load slice profile, this needs to be
                            %%% available, and then is interpolated to
                            %%% different resolution?
                            load(param.slice_profile_file)
                            if ~strcmp(param.resolution,'0p4') 
                                %%% interpolate in space
                                for sim_fa = 1:size(Sprofile,2) %%% for each column
                                    Slab = Sprofile(:,sim_fa);
                                    Slab2(:,sim_fa) = interp1(1:length(Slab),Slab,linspace(1,length(Slab),size(PDw_echo1.data,3)));
                                end
                                Sprofile = Slab2;
                            end
                            [A_echo1_spc,R1_echo1_spc,b1_corr_fa_maps,actual_fa_maps,R1_error] = weighted2AR1_slice_profile_correction({PDw_echo1,T1w_echo1},TR,'spr',percB1/100,Sprofile,90*(0.05:0.05:1));
                            MAPS.R1_echo1_spc = R1_echo1_spc.data;
                            [A_TE0_WLS_spc,R1_TE0_WLS_spc,b1_corr_fa_maps,actual_fa_maps,R1_error] = weighted2AR1_slice_profile_correction({PDw0_WLS,T1w0_WLS},TR,'spr',percB1/100,Sprofile,90*(0.05:0.05:1));
                            MAPS.R1_TE0_WLS_spc = R1_TE0_WLS_spc.data;
                            %%% estimate b1+ for unselective                          
                            pdw_b1rel_total = 1 + (actual_fa_maps(:,:,:,1) - rad2deg(weighted.data.pdw.fa)) ./ rad2deg(weighted.data.pdw.fa);
                            pdw_b1rel_spc = 1 + (pdw_b1rel_total - percB1/100);

                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% calculate MT stuff if MTw exists
                    if isfield(weighted.data, 'mtw') %exist('weighted.data.mtw') ~= 0
                        approaches = {'echo1','TE0_WLS'};
                        for a = 1:length(approaches)
                            approach = approaches{a};
                            if strcmp(approach,'echo1')
                               MTw_to_use = averageEchoes(weighted.data.mtw,1);
                               PDw_to_use = PDw_echo1;
                               if param.mtsatusesmallFA == 1 %%% if the small FA assumption is used
                                R1_to_use = R1_echo1_b1uncorr; %%% conventional MTsat calculation needs b1 uncorrected maps as input
                                A_to_use = A_echo1_b1uncorr;
                               else %%% this is now the new default, to use the corrected maps as input and our newly developed b1 correction for mtsat
                                R1_to_use = R1_echo1; 
                                A_to_use = A_echo1; 
                               end
                            elseif strcmp(approach,'TE0_WLS')
                               MTw_to_use = MTw0_WLS;
                               PDw_to_use = PDw0_WLS;
                               if param.mtsatusesmallFA == 1
                                   R1_to_use = R1_TE0_WLS_b1uncorr;
                                   A_to_use = A_TE0_WLS_b1uncorr;
                               else
                                   R1_to_use = R1_TE0_WLS; 
                                   A_to_use = A_TE0_WLS; 
                               end
                            end
                            %%% MTR
                            mapname = strcat('MTR_',approach);
                            MAPS.(mapname) = (PDw_to_use.data-MTw_to_use.data)./(PDw_to_use.data);
                            MAPS.(mapname)(isnan(MAPS.(mapname))) = 0;
                            %%% MTsat
                            mapname = strcat('MTsat_real_',approach);
                            if param.mtsatusesmallFA == 1 %%% old approach
                                MTsat = computeMTsatApprox(A_to_use,R1_to_use,MTw_to_use,TR);
                            else %%% new approach
                                MTsat = computeMTsatApprox_using_exact_FAs(A_to_use,R1_to_use,MTw_to_use,TR,percB1);
                                %%% if i divide this by ft^2, then i get a
                                %%% map that is "intrinsically corrected"
                            end
                            MAPS.(mapname) = MTsat.data;
                            %%% slice profile correction
                            if (strcmp(param.resolution,'0p4') || strcmp(param.resolution,'0p3')) && param.dosliceprofilecorrection == 1
                                MTsat_spc = pdw_b1rel_spc.^2 .* MAPS.(mapname);
                                mapname = [mapname,'_spc']; 
                                MAPS.(mapname) = MTsat_spc; 
                            end
                            %% rescale to assumed mt pulse flip angle of 700
                            if param.rescale_mt_pulse_degree > 0 %%% this specifies target flip angle
                                first_mtw_filename = [weighted.filenames.mtw(1).folder,'/',weighted.filenames.mtw(1).name];
                                first_mtw_json = [first_mtw_filename(1:end-3),'json'];
                                [desired_fa, sinc_pulse_duration, srf01, srf02, mt_gaussian] = get_fa_parameters_from_json_file(first_mtw_json);
                                nominal_fa_rad = deg2rad(mt_gaussian);
                                desired_fa_rad = deg2rad(param.rescale_mt_pulse_degree);
                                %MTrescaled = interpolate_mt_map_to_different_mt_pulse_with_C(MAPS.(mapname), param.mt_b1_correction_slope, nominal_fa_rad, desired_fa_rad);
                                fake_b1map = (percB1 .^ 0) * 100 * (nominal_fa_rad / desired_fa_rad); %%% fake consistent b1map in perc
                                MTrescaled = b1_correct_mt_map_linear(MAPS.(mapname), fake_b1map, param.mt_b1_correction_C); %%% new correction function
                                mapname1 = [mapname,'_rescaled_to_FA',num2str(param.rescale_mt_pulse_degree)]; 
                                MAPS.(mapname1) = MTrescaled;   
                            else
                                mapname1 = mapname;
                            end
                            %% then b1+ correct rescaled or normal map
                            MTb1corr = b1_correct_mt_map_linear(MAPS.(mapname1), percB1, param.mt_b1_correction_C);
                            mapname1 = strcat(mapname1,'_b1corrected');
                            MAPS.(mapname1) = MTb1corr;
                        end

                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% save and brain mask all in the MAPS structure
                    
                     save_and_brain_mask_files_from_MAPS_structure(MAPS, MPMvol, brain_mask_file, mpm_output_dir, processed_file_string_with_echos, mpm_run, param.resolution)
                     
                     %% this was an attempt to interpolate the values in the buubles with surrounding values, have never optimised
                     if param.createbubblemask == 1
                         %%% this part is WIP for interpolating from
                         %%% surrounding areas
                         contrasts_to_interpolate = {'R2s_OLS','R2s_WLS'};
                         for c = 1:length(contrasts_to_interpolate)
                            mapname = contrasts_to_interpolate{c};
                            brain_masked_filename = [mpm_output_dir, '/', mapname,'_', processed_file_string, '_', param.resolution,'_run', sprintf('%.02d',mpm_run), '_brain_masked.nii'];
                            bubble_mask_preliminary_filename = [mpm_output_dir, '/bubble_mask_preliminary_', processed_file_string, '_', param.resolution,'_run', sprintf('%.02d',mpm_run), '_brain_masked.nii'];
                            bubble_mask_final_filename = [mpm_output_dir, '/bubble_mask_', processed_file_string, '_', param.resolution,'_run', sprintf('%.02d',mpm_run), '_brain_masked.nii'];
                            brain_masked_interpolated_filename = [mpm_output_dir, '/', mapname,'_interpolated_', processed_file_string, '_', param.resolution,'_run', sprintf('%.02d',mpm_run), '_brain_masked.nii'];
                            %%% try to improve the bubble mask a bit
                            system(['bash ~/Documents/scripts/postmortembrain-mpm/Calculation/functions/interpolate_bad_r2s ', brain_masked_filename, ' ', bubble_mask_preliminary_filename, ' ', bubble_mask_final_filename,' ', brain_masked_interpolated_filename, ' ', mpm_output_dir]);
                         end
                         %%% this part fills the bubble gaps eith echo1
                         %%% data
                          contrasts_to_interpolate = {'MTsat_real_TE0_WLS_rescaled_to_FA700_b1corrected','R1_TE0_WLS'};
                          contrasts_to_interpolate_with = {'MTsat_real_echo1_rescaled_to_FA700_b1corrected','R1_echo1'};
                         for c = 1:length(contrasts_to_interpolate)
                            wls_filename = [mpm_output_dir, '/', contrasts_to_interpolate{c},'_', processed_file_string, '_', param.resolution,'_run', sprintf('%.02d',mpm_run), '_brain_masked.nii'];
                            firstecho_filename = [mpm_output_dir, '/', contrasts_to_interpolate_with{c},'_', processed_file_string, '_', param.resolution,'_run', sprintf('%.02d',mpm_run), '_brain_masked.nii'];
                            bubble_mask_final_filename = [mpm_output_dir, '/bubble_mask_', processed_file_string, '_', param.resolution,'_run', sprintf('%.02d',mpm_run), '_brain_masked.nii'];
                            interpolated_filename = [mpm_output_dir, '/', contrasts_to_interpolate{c},'_interpolated_', processed_file_string, '_', param.resolution,'_run', sprintf('%.02d',mpm_run), '_brain_masked.nii'];
                            system(['bash ~/Documents/scripts/postmortembrain-mpm/Calculation/functions/interpolate_WLS_with_first_echo ', wls_filename, ' ', firstecho_filename, ' ', bubble_mask_final_filename,' ', interpolated_filename, ' ', mpm_output_dir]);
                         end
                     end
                    end
                end
            end
        end
    end              
