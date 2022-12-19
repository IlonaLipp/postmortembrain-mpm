param.bashscriptfolder = '/data/hu_lippi/Documents/scripts/postmortembrain-mpm/Calculation/functions' %%% with functions that will be called with this script
param.hmri_processing_parameters =  '/data/hu_lippi/Documents/scripts/postmortembrain-mpm/Calculation/functions/hmri_b1_7T_from_saskia_amended_with_nik_2020_08_25.m'

param.unique_id = 'unknown' %%% sample ID
param.scandate = 000000; %%% in the format YYMMDD, not super relevant for other people than ilona

%% temperature plot
param.dotemp = 0 %%% create temperature plot
param.tempfile = '' %%% filename of file that is stored by truetemp after temperature recording

%% mri files
param.raw_data_folder_nifti = ''; %%% folder where you want your nifti files to end up
param.raw_data_folder_dicom = ''; %%% folder where the dicom files are located
param.redo_nifti_conversion = 1; %%% if you want to do conversion from dicom to nifti, set to 1
param.do_not_check_niftis = 0; %%% if you want to turn nifti checks off (is a feature to redo nifti conversion if something went wrong, does not work with all acquisitions though)
param.do_nifti_conversion_connectom_style = 0; %%% if your dicom directy has structure like from connectom where dicom files are not sorted in subdirectories

%% output folder
param.processed_data_folder = ''; %%% all outputs will be stored there, folder should exist

%% folder_series_numbers: number of scan folder (S*) put in here. if you calculated several runs, the order of runs
%%% if a sequence was not collected, leave empty
param.folder_series_numbers.t1w = []; %%% magnitude folders only, one per repeat of MPMs
param.folder_series_numbers.pdw = []; %%% magnitude folders only, one per repeat of MPMs, vector must match length of t1w 
param.folder_series_numbers.mtw = []; %%% magnitude folders only, one per repeat of MPMs, can be left empty
param.folder_series_numbers.additionalFA = []; %%% magnitude folders of additional flip angles acquired
param.folder_series_numbers.b1_bs = []; %%% folder numbers of bloch-siegert acquisition
param.folder_series_numbers.afi = []; %%% folder numbers of afi acquisitions
param.folder_series_numbers.mafi = []; %%% folder numbers of afi acquisitions
param.folder_series_numbers.al_b1 = []; %%% for SESTE, kp_seste
param.folder_series_numbers.al_b0 = []; %%% for SESTE gre_field_mapping, for each run in b1 folder, it expects two numbers in b0 folder (magnitude and phase)
param.folder_series_numbers.receive_bias_receive = []; %%% magnitude folders of receive bias scans with receive coil
param.folder_series_numbers.receive_bias_transmit = []; %%% magnitude folders of receive bias scans with transmit coil
param.mp2rage_folder_series_numbers = []; %%% not yet relevant
        
%% specify MPM data to analyse
param.resolution = ''; %%% for MPM, will be used in filenames and for resampling of bias maps (e.g. 0.7mm isotropic would be "0p7" in this string field)

%% masks
param.bm_thr = 100; %%% image intensity threshold for creating brain mask (is applied to first echo of t1w)
param.own_brain_mask_file = ''; %%% if you put a nifti filename in here, then it will take this as a brain mask. needs to be in same resolution as weighted images!
%%% and the mask of which will be automatically created

%% B1+ mapping
param.b1_method = {'AL'}; %%% options: 'AL' for Antoine Lutti's SESTE, 'AFI' for Nik's AFI sequence and 'BS' for bloch-siegert, 'custom' for custom b1map that is already calculated (needs to be in percent units!!!) or 'none' to not do b1+ mapping (it will still create a map with all voxels being set to 100); will calculate all you have specified here
param.b1_method_to_use = 'AL' %%% the one to use for MPM calculation, can be set to 'uncorrected' to not do b1+ correction
param.reference_T1 = 500; %%%realistic for postmortem 7T; only relevant for SESTE, assumed T1 in the equation, does not have a strong influence on outcome unless your value is way way lower
param.b1map_smooth_in_mm = 8; %%% smoothing kernel for b1 map
param.custom_b1_filename = ''; %%% if a custom b1 map should be used, provide .nii filename here, and use 'custom' in option b1_method and optionally also b1_method_to_use

%% motion correction
param.do_movement_correction = 0; %%% do only use when sample really fundamentally moved between the weighted sequences
param.mocodof = 12; %%% flirt linear registration is done with the degrees of freedome used here
param.echo2reg = 1; %%% 1 = use first echo for calculating registration matrix

%% distortion correction for distortions between odd and even echos
param.do_distortion_correction = 1; %%% set to 1 if you want distortion correction, it will create a folder processed_s_dc with your distortion corrected data
param.permdim = 1; %%% param for apply_niks_distortion_correction.m

%% receive bias correction
param.do_receive_bias_correction = 0; %%% receive bias corrects all weighted images and saves them in a separate folder
param.rbmap_smooth_in_mm = 8; %%% in Papp 2016 reported as 12 mm, with a 4 mm param.resolution acquisition; in hMRI toolbox it is 7mm kernel? 
param.rbmap_smooth_in_mm_median = 2; %%% additional median filter for smoothing

%% MPMS
param.mtv_calc = 0; %%% set to 0 for no calculation of MTV, and to 1 if you want MTV calculated, this relies on a reference that was scanned in the same session...
param.mtv_refdensity = 0.2000; %%% only relevant if param.mtv_calc is set to 1, density of reference material (1 would be water, if you do not have water, then you need...
%%% to calibrate your reference against water
param.dosliceprofilecorrection = 0; %%% will unlikely be useful, but is still an option
param.slice_profile_file = ''; %%% this needs to be provided if slice profile should be correctd
param.mt_b1_correction_C = 1.2; %%% based on Lipp et al., 2022, MRM
param.rescale_mt_pulse_degree = 700; %%% flip angle to which MTsat maps should be calibrated to, is chimp default
param.echos_to_use = {'all'}; %%% number of echos to use for calculations (min = 6), can have several entries as will loop around, if set to "skipfirst", it gets rid of first echo (in case of saturation)
param.createbubblemask = 0; %%% will try to create a mask for bubbles, will slow down the script massively
param.correct_fa = 1; %%% if set to yes, SRF is checked in header and in case of clipped pulse, actual flip angle is calculated and used
param.mtsatusesmallFA = 0; %%% keep to 0 if you want mtsat maps as described in lipp et al. 2022 MRM, set to 1 if you have low flip angles, and it will use approximation

%% MP2RAGE %%% probably not working
param.processmp2rage = 0; 
param.qmrlabpath = '/data/tu_lippi/Software/qMRLab';
param.Guiopts = [param.bashscriptfolder,'/mp2rage.qmrlab.Magnetom.qmrlab.mat']; 
param.bm_thr_mp2rage = 150; %%% threshold for UNI2 image to get brain mask
param.folder_series_numbers.mp2rage_INV2 = []; %% 1 
param.folder_series_numbers.mp2rage_UNI = []; %% 1 
param.folder_series_numbers.mp2rage_T1 = []; %% 1
param.useforT1ref = 'qmrlab' %%% options are scanner or qmrlab or 'reference