function apply_niks_distortion_correction(input_echos, output_geom, output_bfield, output_corr, permdim)
    %%% writes SPM batches to run distortion correction as implemented by
    %%% nik. it puts the corrected files in a now folder and also copies
    %%% the original json files over, so that they can be loaded into MPM
    %%% processing pipeline as usual
    
%     %%% make temporary dir
%     temp_dir = [output_corr,'/temp'];
%     system(['mkdir ', temp_dir);
    
    spm_jobman('initcfg');
    spm('defaults', 'fMRI');

    imcalmatlabbatch{1}.spm.util.imcalc.input=input_echos([1 3]); %' {input_echos{1};input_echos{3}};
    imcalmatlabbatch{1}.spm.util.imcalc.output=output_geom;
    imcalmatlabbatch{1}.spm.util.imcalc.outdir={output_corr};
    imcalmatlabbatch{1}.spm.util.imcalc.expression='exp((log(i1)+log(i2))/2)';
    imcalmatlabbatch{1}.spm.util.imcalc.var=[];
    imcalmatlabbatch{1}.spm.util.imcalc.options.dmtx=0;
    imcalmatlabbatch{1}.spm.util.imcalc.options.mask=0;
    imcalmatlabbatch{1}.spm.util.imcalc.options.interp=-7;
    imcalmatlabbatch{1}.spm.util.imcalc.options.dtype=4;

    spm_jobman('run',imcalmatlabbatch);
    
    %%% http://www.diffusiontools.com/documentation/hysco.html
    %  HySCO2_main(PI1,PI2,POI1,POI2,pe_direction,full_res,doECC,alpha,beta,dummy_3dor4d)
    %  PI1          - filename of reference blip-up
    %  PI2          - filename of reference blip-down
    %  POI1         - matrix of filenames for additional blip-up volumes
    %  POI2         - matrix of filenames for additional blip-down volumes
    %  pe_direction - phase-encoding direction, 1 for x_1, 2 for x_2, 3 for x_3
    %                 (data dimensions will be flipped accordingly)
    %  full_res     - finest level for multi-level correction, boolean
    %  doECC        - do nonlinear eddy-correction for additinal volumes (equal
    %                 number of blip-up and blip-down data required.)
    %  alpha        - regularization paramter for diffusion term
    %  beta         - regularization paramter for Jacobian term
    %  dummy_3dor4d - Specifies output format
    
    spm_jobman('initcfg');
    spm('defaults', 'fMRI');
      
    matlabbatch{1}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.source_up = {[input_echos{2},',1']}; %{'/data/pt_02101/raw/007_C_C_NEGRA_ID/mr/191026_Magnetom_7T_32Ch_WB/nifti/pdw_kp_mtflash3d_v1d_2p1_0023/s2019-10-26_09-32-104717-00001-00049-2.nii,1'};
    matlabbatch{1}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.source_dw = {[output_geom,',1']}; %{'/data/pt_02101/preprocessed/007_C_C_NEGRA_ID/mr/191026_Magnetom_7T_32Ch_WB/distortion_correction/pdw_kp_mtflash3d_v1d_2p1_0023_geom.nii,1'};
    matlabbatch{1}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.others_up = {''};
    matlabbatch{1}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.others_dw = {''};
    matlabbatch{1}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.perm_dim = permdim;
    matlabbatch{1}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.dummy_fast = 1;
    matlabbatch{1}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.dummy_ecc = 0;
    matlabbatch{1}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.alpha = 50;
    matlabbatch{1}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.beta = 10;
    matlabbatch{1}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.dummy_3dor4d = 0;
    matlabbatch{1}.spm.tools.dti.prepro_choice.hysco_choice.hysco2.restrictdim = [1
                                                                                  1
                                                                                  1];
    % function HySCO_write(POI1,POI2,Bc,pe_direction)
    %  POI1         - matrix of filenames for additional blip-up volumes
    %  POI2         - matrix of filenames for additional blip-down volumes
    %  VB           - filename of inhomogeneity estimate produced by HySCO
    %  pe_direction - phase-encoding direction, 1 for x_1, 2 for x_2, 3 for x_3
    %  (data dimensions will be flipped accordingly)
    %%% add volume number to input echos:
    count = 1;
    for i = 2:2:length(input_echos)
        input_array_even{count} = [input_echos{i},',1'];
        count = count + 1;
    end
    count = 1;
    for i = 1:2:length(input_echos)
        input_array_odd{count} = [input_echos{i},',1'];
        count = count + 1;
    end
    matlabbatch{2}.spm.tools.dti.prepro_choice.hysco_choice.hysco_write.source_up = {[input_echos{2},',1']};
    matlabbatch{2}.spm.tools.dti.prepro_choice.hysco_choice.hysco_write.others_up = input_array_even';
    matlabbatch{2}.spm.tools.dti.prepro_choice.hysco_choice.hysco_write.others_dw = input_array_odd';
    matlabbatch{2}.spm.tools.dti.prepro_choice.hysco_choice.hysco_write.bfield = {[output_bfield,',1']};
    matlabbatch{2}.spm.tools.dti.prepro_choice.hysco_choice.hysco_write.perm_dim = permdim;
    matlabbatch{2}.spm.tools.dti.prepro_choice.hysco_choice.hysco_write.dummy_3dor4d = 0;
   
    [fol fi] = fileparts(input_echos{1});
    system(['yes | rm ', fol, '/u2*2.*']) %%% remove the weird output
    system(['yes | rm ', output_corr, '/u2*2.*']) %%% remove the weird output
    
    %if length(dir(output_corr)) < 5
    spm_jobman('run',matlabbatch);
    %end
    
    system(['yes | rm ', fol, '/u2*2.*']) %%% remove the weird output
    system(['yes | rm ', output_corr, '/u2*2.*']) %%% remove the weird output
    
    %%% move output files
 
    system(['mv ', fol, '/u*.nii ',output_corr])


    %input_echos = convert_dir_output_to_cell_structure(dir([fol,'/*.nii']));
    
    all_outputs = dir([output_corr,'/u*.nii']);
    if length(input_echos) ~= length(all_outputs)
        display('proooblem')
    end
    
    %%% cp json too
    for i = 1:length(input_echos)
        %%% get name of json file
        [fol fi] = fileparts(input_echos{i});
        old_json_file = [fol,'/',fi,'.json'];
        %%% get name of output image
        [fol fi] = fileparts([all_outputs(i).folder,'/',all_outputs(i).name]);
        new_json_file = [fol,'/',fi,'.json'];
        %%% are same alphabetical order
        system(['cp ',old_json_file,' ',new_json_file]);
    end
end