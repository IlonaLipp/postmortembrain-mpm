function seriesinfo = analyse_scan_parameters_from_json_files_for_a_scan_session(param,session_folder, pdfoutputfilename, frequencyplot_outputfilename)
        %%% this function goes through all folders in a scan session (nifti
        %%% folders as converted (with hmri toolbox) from dicom folder
        %%% saved by Siemens) and reads information in the json files to
        %%% make a pdf summary of how things were acquired
        %%% session_folder = folder_name
        %%% pdfoutputfilename = file where information should be written
        %%% frequencyplot_outputfilename = file to write frequency plot
        %%% output is seriesinfo structure that saves some information for
        %%% further processing

        %seriesinfo = 0;
        
        if exist(pdfoutputfilename) ~= 0
        delete(pdfoutputfilename); %%% because written with -append
        end
        scanfoldernames = dir(session_folder);
        if length(scanfoldernames) > 0
            clear seq_number
            %%% sort by when scanned
            for s = 1:length(scanfoldernames)
                scanfoldername = scanfoldernames(s).name;
                split = strsplit(scanfoldername,'_');
                seq_number(s) = str2double(split{end});
            end
            [sorted idx] = sort(seq_number);
            scanfoldernames = scanfoldernames(idx);
            clear frequency_drift
            frequency_drift_time = 0;
            scanfoldernames
            for s = 1:length(scanfoldernames) %%% first two are just ., and ..

                sequence_number = sorted(s)
                scanfoldername = [scanfoldernames(s).folder, '/', scanfoldernames(s).name];
                filenames = dir([scanfoldername,'/*.json']);
                tableinfo = {};

                %%% only print certain
                if length(filenames) > 0 && count(scanfoldernames(s).name, '_') > 0 && count(scanfoldername, 'RR') == 0 %%&& count(scanfoldername, 'mtflash3d') > 0 && (count(scanfoldername, '0p3') > 0 || count(scanfoldername, '0p4') > 0  || count(scanfoldername, '400um') > 0 || count(scanfoldername, '1p6') > 0) %|| count(scanfoldername, 'mtflash3d')
                    display('passed first test')
                    clear TEs
                    for f = 1:length(filenames) %%% loop over echos just to get all echo info
                        filename = filenames(f).name;
                        json_chars = fileread([scanfoldername,'/',filename]); %%% read associated json file 
                        json_chars = strrep(json_chars,',]',']'); %%% jsondecode struggles with ',]', so remove comma in all occurences
                        json = jsondecode(json_chars); %%% read associated json file 
                        %%% read first echo
                        TEs(f) = json.acqpar.EchoTime; %%% save echo times in ms in vector
                        %%% use information from first echo as start of scan
                        if f == 1 %%% base date and time on first echo in each sequence
                            seriesinfo(sequence_number).name = json.acqpar.SeriesDescription;
                            date = datestr(json.acqpar.AcquisitionDate); %%% is stored as days from 0
                            hours_first = floor(json.acqpar.AcquisitionTime/3600);
                            minutes_first = floor((json.acqpar.AcquisitionTime/3600 - hours_first) * 60);
                            seconds_first = floor(60 * ((json.acqpar.AcquisitionTime/3600 - hours_first) * 60 - minutes_first));
                            %%% store in seriesinfo
                            try
                                seriesinfo(sequence_number).time = [sprintf('%02d',hours_first),':',sprintf('%02d',minutes_first),':',sprintf('%02d',seconds_first)];
                            end
                        end
                        if f == length(filenames)
                            hours_last = floor(json.acqpar.AcquisitionTime/3600);
                            minutes_last = floor((json.acqpar.AcquisitionTime/3600 - hours_last) * 60);
                            seconds_last = floor(60 * ((json.acqpar.AcquisitionTime/3600 - hours_last) * 60 - minutes_last));
                        end
                    end         

                    %%% only proceed if the same time has not been used just
                    %%% now
                    display('displaying filename')
                    seriesinfo.name
                  %  if s > 1 && length(seriesinfo(s-1).name) > 6 && (isequal(seriesinfo(s).name(1:7),seriesinfo(s-1).name(1:7)) && isequal(seriesinfo(s).time, seriesinfo(s-1).time))
                    if s > 1 && isequal(seriesinfo(sorted(s)).time, seriesinfo(sorted(s)-1).time)
                    %if sequence_number > 1 && isequal(seriesinfo(s).time, seriesinfo(sequence_number-1).time)    
                        display('same time, not writing out anything')
                        seriesinfo(sequence_number).time;
                        seriesinfo(sequence_number-1).time;
                    else
                         display('writingout')
                        %%% use information from last echo to write out anything else
                        tableinfo = cat(1, tableinfo, {'Series number', num2str(sequence_number)});
                        tableinfo = cat(1, tableinfo, {'Foldername',scanfoldernames(s).name});
                        tableinfo = cat(1, tableinfo, {'AcquisitionDate',date}); 
                      %  tableinfo = cat(1, tableinfo, {'AcquisitionTime',[sprintf('%02d',hours_first),':',sprintf('%02d',minutes_first)]}); %%% is sec from 0:00
                       % tableinfo = cat(1, tableinfo, {'AcquisitionTime last',[sprintf('%02d',hours_last),':',sprintf('%02d',minutes_last)]}); %%% is sec from 0:00
                        tableinfo = cat(1, tableinfo, {'AcquisitionTime',seriesinfo(s).time});
                        tableinfo = cat(1, tableinfo, {'SeriesDescription', json.acqpar.SeriesDescription}); 
                        tableinfo = cat(1, tableinfo, {'StudyDescription', json.acqpar.StudyDescription});
                        tableinfo = cat(1, tableinfo, {'SequenceName', json.acqpar.SequenceName});
                        tableinfo = cat(1, tableinfo, {'MagneticFieldStrength (T)', num2str(json.acqpar.MagneticFieldStrength)});
                        tableinfo = cat(1, tableinfo, {'ImagingFrequency (Hz)', num2str(json.acqpar.ImagingFrequency)});
                        tableinfo = cat(1, tableinfo, {'NumberOfPhaseEncodingSteps', num2str(json.acqpar.NumberOfPhaseEncodingSteps)});
                        tableinfo = cat(1, tableinfo, {'PercentSampling', num2str(json.acqpar.PercentSampling)});
                        tableinfo = cat(1, tableinfo, {'PercentPhaseFieldOfView', num2str(json.acqpar.PercentPhaseFieldOfView)});
                        tableinfo = cat(1, tableinfo, {'PixelBandwidth (Hz)', num2str(json.acqpar.PixelBandwidth)});
                        tableinfo = cat(1, tableinfo, {'SliceThickness (mm)', num2str(json.acqpar.SliceThickness)});
                        tableinfo = cat(1, tableinfo, {'In plane resolution (mm)', [num2str(json.acqpar.PixelSpacing(1)),' x ',num2str(json.acqpar.PixelSpacing(2))]});
                        tableinfo = cat(1, tableinfo, {'AcquisitionMatrix size', num2str(json.acqpar.AcquisitionMatrix')});
                        tableinfo = cat(1, tableinfo, {'RepetitionTime (ms)',  num2str(json.acqpar.RepetitionTime)});
                        tableinfo = cat(1, tableinfo, {'FlipAngle (deg)', num2str(json.acqpar.FlipAngle)});
                        try
                        tableinfo = cat(1, tableinfo, {'ucPhasePartialFourier', num2str(json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sKSpace.ucPhasePartialFourier)});
                        tableinfo = cat(1, tableinfo, {'ucSlicePartialFourier', num2str(json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sKSpace.ucSlicePartialFourier)});
                        end
                        try
                        tableinfo = cat(1, tableinfo, {'Receive coil', json.acqpar.CSAImageHeaderInfo.ImaCoilString});
                        end
                        try
                        tableinfo = cat(1, tableinfo, {'sPat.lAccelFact3D', num2str(json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sPat.lAccelFact3D)});
                        tableinfo = cat(1, tableinfo, {'sPat.lAccelFactPE', num2str(json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sPat.lAccelFactPE)});
                        tableinfo = cat(1, tableinfo, {'FirstEcho (ms)', num2str(TEs(1))});
                        tableinfo = cat(1, tableinfo, {'LastEcho (ms)', num2str(TEs(end))});
                        tableinfo = cat(1, tableinfo, {'#Echos', num2str(length(TEs))});
                        tableinfo = cat(1, tableinfo, {'TransmitVoltage', num2str(json.acqpar.CSASeriesHeaderInfo.TransmitterCalibration)});
                        end
                        %%% MPM specific
                        try %%% only works if the field exists
                            try
                                sinc_pulse_duration = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sWiPMemBlock.alFree(8);
                            catch %%% new scanner
                              sinc_pulse_duration = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sWipMemBlock.alFree(17);
                              json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sWipMemBlock.alFree';
                              json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sWipMemBlock.adFree';
                            end
                            tableinfo = cat(1, tableinfo, {'sinc pulse dur (us)', num2str(sinc_pulse_duration)});  
                        end
                        try
                            try %%% old scanner
                                mt_gaussian = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sWiPMemBlock.adFree(1);
                                mt_pulse_duration = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sWiPMemBlock.alFree(1);
                            catch %%% new scanner
                                mt_gaussian = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sWipMemBlock.adFree(1);
                                mt_pulse_duration = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sWipMemBlock.alFree(2)
                            end
                            try
                                mt_volt = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sTXSPEC.aRFPULSE(1).flAmplitudeNL; %% new scanner
                                srf01 = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sTXSPEC.aRFPULSE(2).flAmplitudeNL;
                                srf02 = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sTXSPEC.aRFPULSE(3).flAmplitudeNL;
                            catch
                               mt_volt = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sTXSPEC.aRFPULSE(1).flAmplitude;
                                srf01 = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sTXSPEC.aRFPULSE(2).flAmplitude;
                                srf02 = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sTXSPEC.aRFPULSE(3).flAmplitude;   
                            end
                            %%% write in table
                            tableinfo = cat(1, tableinfo, {'sSRF01 (V)', num2str(srf01)});        
                            tableinfo = cat(1, tableinfo, {'sSRF02 (V)', num2str(srf02)}); 
                            if isfield(json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sPrepPulses,'ucMTC')
                                tableinfo = cat(1, tableinfo, {'MT', 'on'});
                                tableinfo = cat(1, tableinfo, {'MTGaussian', num2str(mt_gaussian)}); 
                                tableinfo = cat(1, tableinfo, {'MTVolt', num2str(mt_volt)}); 
                                tableinfo = cat(1, tableinfo, {'MTPulseDuration', num2str(mt_pulse_duration)}); 
                            else
                                tableinfo = cat(1, tableinfo, {'MT', 'off'});
                            end
                        end
                        %%% shim pars
                        try
                        if json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sAdjData.uiAdjShimMode == 1
                            tableinfo = cat(1, tableinfo, {'Shim type', 'tuneup'});
                        elseif json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sAdjData.uiAdjShimMode == 2
                            tableinfo = cat(1, tableinfo, {'Shim type', 'not tuneup'});
                        end
                        tableinfo = cat(1, tableinfo, {'Shim pars', num2str(json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sGRADSPEC.alShimCurrent')});
                        end

                        %%% coil combination
                        try
                        if json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.ucCoilCombineMode == 2
                            tableinfo = cat(1, tableinfo, {'Coil combination', 'adaptive combine'});
                        elseif json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.ucCoilCombineMode == 1
                            tableinfo = cat(1, tableinfo, {'Coil combination', 'SSQ'});
                        end
                        end

                        %%% print for now
                        json.acqpar.SeriesDescription
                        if json.acqpar.AcquisitionDate * 24 * 60 * 60 + json.acqpar.AcquisitionTime > frequency_drift_time(end) %% just unique entries
                            frequency_drift(s) = json.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sTXSPEC.asNucleusInfo(1).lFrequency;
                            frequency_drift_time(s) = json.acqpar.AcquisitionDate * 24 * 60 * 60 + json.acqpar.AcquisitionTime;
                            frequency_series_name{s} = json.acqpar.SeriesDescription;
                        end

                        %% get maximum intensity
                        [error,max_int] = system(['bash ',param.bashscriptfolder,'/check_sequence_for_saturation ',scanfoldername]);
                        max_int = str2double(max_int);
                        tableinfo = cat(1, tableinfo, {'95th percentile image intensity', num2str(max_int)});

                        %% mt pulse, matrix, bandwidth, partial fourier, grappa, ramp, offset
                        close all
                        f = figure(1);
                        set(gcf,'Visible','off');
                        set(f,'position',[.3 .3 .4 .6],'units','normalized'); %,[100 100 800 1200])
                        set(f,'position',[.3 .1 .4 .8],'units','normalized'); 
                        uit = uitable(f); %, 'Data', tableinfo, 'Position', [0 0 1 1]);
                        %uit.ColumnWidth = {200,500};
                        uit.Units = 'normalized';
                        res = get(0,'screensize');
                        uit.ColumnWidth = {res(3)*.15,res(3)*.25}; % {'auto','auto'} %;{50,50};
                        uit.Position = [0 0 1 1] ;%[0 0 700 1100];
                        uit.Data = tableinfo;
                        %set(f,'position',[.3 .3 .4 .6],'units','normalized');
                        set(f,'position',[.3 .1 .4 .8],'units','normalized');
                        print(f,'-dpsc','-append',pdfoutputfilename);
                    end
                %else
                    %seriesinfo(sequence_number).name = 'unsure';
                    %seriesinfo(sequence_number).time = [];
                %end
                end
            end

        if exist('frequency_drift') ~= 0
            idx = frequency_drift>0;
            frequency_drift_time = frequency_drift_time - frequency_drift_time(1);
            frequency_drift_delta = frequency_drift - frequency_drift(1);
            f = figure(2);
            set(gcf,'Visible','off');
            set(f,'position',[0 0 2000 1000]);
            plot(frequency_drift_time(idx), frequency_drift_delta(idx),'-k');
            xticks((frequency_drift_time(idx)));
            xticklabels((frequency_series_name(idx)));
            set(gca,'TickLabelInterpreter','none');
            xtickangle(80);
            xlabel('scan');
            ylabel('Hz change from time 0');
            saveas(f, frequencyplot_outputfilename);
            %min(frequency_drift(frequency_drift>0))
            %max(frequency_drift(frequency_drift>0))
        close all
        end
     end
end
