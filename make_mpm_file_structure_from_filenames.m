function mpm_file_structure = make_mpm_file_structure_from_filenames(filenames)
%%% this function gets a list of nifti filenames and makes a matlab
%%% structure out of there content for Luke's function to run


if ~isempty(filenames)
    VCON = spm_vol(fullfile(filenames(1).folder,filenames(1).name)); %%% read in first file to define dims
    CON.data = NaN([VCON.dim,numel(filenames)]); %%% initiate volume data; changed to NaN from zeros
    CON.TEs = NaN([1,numel(filenames)]); %%% initiate TEs
    for n = 1:numel(filenames)
        Vtmp = spm_vol(fullfile(filenames(n).folder,filenames(n).name)); %%% read volume of file in loop
        %filenames(n).name
        size(spm_read_vols(Vtmp));
        CON.data(:,:,:,n) = spm_read_vols(Vtmp); %%% concatenate all volumes
        [prefix,name] = fileparts(filenames(n).name); %%% remove extension
        json_chars = fileread(fullfile(filenames(n).folder,prefix,[name,'.json'])); %%% read associated json file 
        json_chars = strrep(json_chars,',]',']'); %%% jsondecode struggles with ',]', so remove comma in all occurences
        json = jsondecode(json_chars); %%% read associated json file 
        CON.TEs(n) = json.acqpar.EchoTime; %%% save echo times in ms in vector
    end
    CON.fa = deg2rad(json.acqpar.FlipAngle); %%% get flip angle and convert to radians, which is what lukes calculations need
    CON.TR = json.acqpar.RepetitionTime*1e-3; %%% get TR in s
    %CON = selectEchoes(CON,1:json.acqpar.EchoNumbers); %%% from concatenated volumes remove phase data
else 
    CON = 0;
end
mpm_file_structure = CON;