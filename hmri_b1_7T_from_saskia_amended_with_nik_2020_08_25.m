function hmri_b1_7T
% Sets the defaults for B1 bias correction, part of the hMRI toolbox.
% Consider this file as a template for local settings specifications. 
% Please read below for details.
%
% FORMAT hmri_b1_local_defaults
%__________________________________________________________________________
%
% PURPOSE
% To set user-defined (site- or protocol-specific) defaults parameters for
% B1 mapping. Applies to 3D EPI, 3D AFI and UNICORT protocols only. 
% Customized processing parameters can be defined, overwriting defaults
% from hmri_b1_standard_defaults. Acquisition parameters can be specified
% here as a fallback solution when no metadata are available. Note that the
% use of metadata is strongly recommended.  
%
% RECOMMENDATIONS
% Parameters defined in this file are identical, initially, to the ones
% defined in hMRI-Toolbox\config\hmri_b1_standard_defaults. It is
% recommended, when modifying this file, to remove all unchanged entries
% and SAVE THE MODIFIED FILE WITH A MEANINGFUL NAME. This will help you
% identifying the appropriate defaults to be used during map creation for
% B1 map calculation, and will improve the readability of the file by
% pointing to the modified parameters only. 
%
% WARNING
% Modification of the defaults parameters may impair the integrity of the
% toolbox, leading to unexpected behaviour. Only recommended for expert
% users. 
%
% HOW DOES IT WORK?
% The modified defaults file can be selected when specifying the B1 type in
% the "Create maps" branch of the hMRI-Toolbox.
%__________________________________________________________________________
% Written by E. Balteau, 2017.
% Cyclotron Research Centre, University of Liege, Belgium
%__________________________________________________________________________

% Global hmri_def variable used across the whole toolbox
global hmri_def

%--------------------------------------------------------------------------
% B1 mapping processing parameters 
%--------------------------------------------------------------------------
% Default parameters are set below for each type of B1 processing.
% For acquisition parameters, default values are a fallback solution for B1
% data processing when no metadata are available. Use of metadata is
% recommended to retrieve site- & protocol-specific parameters and ensure
% appropriate data handling and processing.
% See examples of local customization in the hMRI-Toolbox\local directory.
% 'i3D_EPI'
% hmri_def.b1map.i3D_EPI.b1type = '3D_EPI_v2f_long_7T34910'; 
hmri_def.b1map.i3D_EPI.b1type = 'i3D_EPI'; 
% hmri_def.b1map.i3D_EPI.deffnam  = {fullfile(fileparts(mfilename('fullpath')),'config','local','hmri_b1_local_defaults.m')};
% b0&b1-processing
hmri_def.b1map.i3D_EPI.b1proc.T1 = 1633; % in ms, strictly valid only at 3T, assumed T1 for the b1 calculation, does not have a strong effect unless it is very far off
hmri_def.b1map.i3D_EPI.b1proc.HZTHRESH = 190; %%% these are for processing, ilona's script does it's own processing
hmri_def.b1map.i3D_EPI.b1proc.SDTHRESH = 6; %%% these are for processing, ilona's script does it's own processing
hmri_def.b1map.i3D_EPI.b1proc.B1FWHM = 8;  %%% these are for processing, ilona's script does it's own processing % For smoothing. FWHM in mm - i.e. it is divided by voxel resolution to get FWHM in voxels
% b1-acquisition: these are parameters that are dependent on the specific sequence implementation and need to be checked with whoever set up the sequence
hmri_def.b1map.i3D_EPI.b1acq.beta = 165:-7.5:60; %%% So, the sequence special card / json file might say 330 in 15 steps of 15 degrees. In the toolbox this is entered as 165 in 15 steps of 7.5 degrees. (165:-7.5:60?) What beta actually is: the sequence is repeated 15 times and each time a spin echo ("beta" degrees RF pulse) and stimulated echo (beta/2 degrees RF pulse) are acquired. There are 15 sets of beta so that we can be sure we can reliably calculate B1 from the ratio between SE & STE across the field of view, no matter whether the B1 is high or low. The toolbox processing works out which betas to trust at each point in space.

hmri_def.b1map.i3D_EPI.b1acq.TM = 34.91;
% b0-acquisition: these are parameters that are dependent on the specific sequence implementation and need to be checked with whoever set up the sequence
hmri_def.b1map.i3D_EPI.b0acq.longTE = 11.02; % ms
hmri_def.coreg2PDw = 0;


hmri_def.errormaps  = true;
hmri_def.hom        = false;
hmri_def.qMRI_maps_thresh.R2sHO       = 1;    % [1/s^2]
hmri_def.wols       = true;
% these belong to the weighted least square fit
hmri_def.wolsdef.thr_w0 = 0.01; % regularization factor
hmri_def.wolsdef.sigmaMPM = log(50); % noise - should go in via gui
% hmri_def.wolsdef.brainmask =
% %nw: no brainmask
hmri_def.wolsdef.brainmask = [];
end
