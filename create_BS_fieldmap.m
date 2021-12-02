function create_BS_fieldmap(B1files, outputfilenameB1, outputfilenameDiff, maskoutputfilename, smth)
    %%% this function calculates B1+ map BSSB style from acquired nifti
    %%% files (magnitude and phase difference) and saves it and saves a
    %%% mask too. smoothing (smth = kernel) is optional, would not recommend, but instead
    %%% recommend doing home made b1+ processing pipeline
    %%% functino is based on Luke's

    %%% make proforma brainmask
    system(['fslmaths ',B1files(1).folder,'/',B1files(1).name,' -thr 100 -bin ',maskoutputfilename]);

    [prefix,name] = fileparts(B1files(3).name); %%% the third file of this sequence is of interest?
    json_filename = fullfile(B1files(3).folder,prefix,[name,'.json']);
    %json_chars = fileread(fullfile(B1files(3).folder,prefix,[name,'.json'])); %%% this is the associated jsfon file
    %json_chars = strrep(json_chars,',]',']'); %%% jsondecode struggles with ',]', so remove comma in all occurences
    %json = jsondecode(json_chars); %%% read 

    VB1phase = spm_vol(fullfile(B1files(3).folder,B1files(3).name)); %%% third and fourth files alphabetically is always the phase?
    phaseImg1 = spm_read_vols(VB1phase)/4096*pi; %%% read phase image and normalise?
    phaseImg2 = spm_read_vols(spm_vol(fullfile(B1files(4).folder,B1files(4).name)))/4096*pi; %%% read other phase
    BSintensity = spm_read_vols(spm_vol(fullfile(B1files(1).folder,B1files(1).name))); %%% get intensity map
    BSmask = BSintensity > 0; %128; %%% create rough brain mask

    %Vref = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sTXSPEC.asNucleusInfo(1).flReferenceAmplitude; %%% read reference amplitude from header (is this from kerrin's special sequence?)
    %Vpulse = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sTXSPEC.aRFPULSE(1).flAmplitude; %%% read Amplitude
    %gamma = 2.675222005e8; % rad / s / T; 1H gyromagnetic ratio
    %B1ref = pi/(gamma*1e-3); % T; field for pi pulse over 1 ms (Siemens definition)
    %KBS = 4.93661e9;

    %% added 27/03/19 Luke's new function
    [Vref,Vpulse,KBS,B1ref] = readBSSparamsFromJSON(json_filename);
    
    phaseImg1 = phaseImg1 .* BSmask; %%% mask phase images
    phaseImg2 = phaseImg2 .* BSmask;

    B1 = BlochSiegertB1Map(phaseImg2,phaseImg1,KBS,Vpulse,Vref,B1ref); %%% calculate B1

    B1(B1<0|B1>2) = 0; %%% get rid of unrealistic values?
    B1(abs(imag(B1))>0) = 0; %%% ?

    %%% save difference image
    
    %%% save output in processed files folder and smooth
    VB1 = VB1phase; 
    VB1.dt(1) = spm_type('float64');
    VB1.fname = outputfilenameB1;
    spm_write_vol(VB1,real(B1 * 100)); %%% so that units are consistent
    
    phaseDiff=bsxfun(@minus,phaseImg2,phaseImg1);
    VB1.fname = outputfilenameDiff;
    spm_write_vol(VB1,phaseDiff);
    
    %pxs = sqrt(sum(VB1.mat(1:3,1:3).^2)); % Voxel resolution
    %smth = 8./pxs;
    if smth > 0
    spm_smooth(VB1.fname,VB1.fname, smth); %%% smooth with smooth_factor mm FWHM; make comparable to hmri toolbox, which smoothes with pxs = sqrt(sum(V1.mat(1:3,1:3).^2)); % Voxel resolution smth = 8./pxs;
    end
end