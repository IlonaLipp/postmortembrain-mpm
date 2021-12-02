function [T, TR_s, EchoSpace_ms, FA_deg, TI_s, preshots, postshots] = get_mp2rage_parameters_from_json_file(filename)
    
    json_chars = fileread(filename); %%% read associated json file 
    json_chars = strrep(json_chars,',]',']'); %%% jsondecode struggles with ',]', so remove comma in all occurences
    json = jsondecode(json_chars); %%% read associated json file 

    T = json.acqpar.MagneticFieldStrength;
    TR_s = json.acqpar.RepetitionTime / 1000;
    %echosp
    EchoTrainDuration_ms = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sFastImaging.lEchoTrainDuration;
    EchoPartition = json.acqpar.CSAImageHeaderInfo.EchoPartitionPosition; %%% does not make sense but is only value that would get us to 6ms
    EchoTime_ms = json.acqpar.EchoTime;
    EchoSpace_ms = (EchoTrainDuration_ms / (EchoPartition)) / EchoTime_ms;
    
    FA_deg = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.adFlipAngleDegree;
    TI_s = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.alTI / 1000000;
    %%% calculate shots pre post
    slices = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sKSpace.lImagesPerSlab;
    partial_fourier_code = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sKSpace.ucSlicePartialFourier;
    %%% these conversions luke sent me via email, and https://github.com/ismrmrd/siemens_to_ismrmrd/blob/master/parameter_maps/IsmrmrdParameterMap.xsl#L37
    if partial_fourier_code == 4
        partial_fourier = 6/8;
    elseif partial_fourier_code == 1
        partial_fourier = 4/8;
    elseif partial_fourier_code == 2
        partial_fourier = 5/8;
    elseif partial_fourier_code == 8 
        partial_fourier = 7/8; 
    elseif partial_fourier_code == 16
        partial_fourier = 8/8;
    end
    preshots = slices * (partial_fourier - 0.5);
    postshots = slices * 0.5;
end
