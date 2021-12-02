function [desired_fa, tr_diff] = get_mafi_parameters_from_json_file(filename)

    %%% this function reads out information related to the flip angle from
    %%% a json file
    
    json_chars = fileread(filename); %%% read associated json file 
    json_chars = strrep(json_chars,',]',']'); %%% jsondecode struggles with ',]', so remove comma in all occurences
    json = jsondecode(json_chars); %%% read associated json file 


    desired_fa = json.acqpar.FlipAngle;
    tr_diff = json.acqpar.CSAImageHeaderInfo.EchoPartitionPosition; %%% not sure that this is correct!!!
    
end
         
