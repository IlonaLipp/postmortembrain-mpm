function [desired_fa, sinc_pulse_duration, srf01, srf02, mt_gaussian] = get_fa_parameters_from_json_file(filename)

    %%% this function reads out information related to the flip angle from
    %%% a json file
    
    json_chars = fileread(filename); %%% read associated json file 
    json_chars = strrep(json_chars,',]',']'); %%% jsondecode struggles with ',]', so remove comma in all occurences
    json = jsondecode(json_chars); %%% read associated json file 


    desired_fa = json.acqpar.FlipAngle;
    try
        try %%% old scanner:
            sinc_pulse_duration = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sWiPMemBlock.alFree(8);
            mt_gaussian = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sWiPMemBlock.adFree(1);
        catch %%%then it is probably acquired on the new scanner
            sinc_pulse_duration = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sWipMemBlock.alFree(17);
            mt_gaussian = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sWipMemBlock.adFree(1);
        end
        try
            srf01 = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sTXSPEC.aRFPULSE(2).flAmplitudeNL; %% new scanner
            srf02 = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sTXSPEC.aRFPULSE(3).flAmplitudeNL;
        catch
            srf01 = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sTXSPEC.aRFPULSE(2).flAmplitude;
            srf02 = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sTXSPEC.aRFPULSE(3).flAmplitude;
        end
        try %%% old scanner
            mt_gaussian = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sWiPMemBlock.adFree(1);
            mt_pulse_duration = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sWiPMemBlock.alFree(1);
        catch %%% new scanner
            mt_gaussian = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sWipMemBlock.adFree(1);
            mt_pulse_duration = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sWipMemBlock.alFree(2);
        end
    catch
            sinc_pulse_duration = NaN;
            mt_gaussian = NaN;
            srf01 = NaN;
            srf02 = NaN;
            mt_gaussian = NaN;
    end
