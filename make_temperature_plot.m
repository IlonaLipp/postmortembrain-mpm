function make_temperature_plot(filename,seriesinfo,outputfilename_temp, outputfilename_fig, scanner)
        %%% this function makes a plot showing the temperature trajectory
        %%% across a scan session. for this it needs a file that has the
        %%% temperature recordings and time stamp and a seriesinfo
        %%% structure that was created by the function
        %%% analyse_scan_parameters_from_json_files_for_a_scan_session and
        %%% has information about what sequences were run at what time
        %%% filename = filename of the temperature recording
        %%% seriesinfo = structure with series information
        %%% outputfilename_temp = stores temperature information as a
        %%% matlab variable for further analysis
        %%% outputfilename_fig = filename to store figure
        %%% DIFFERENCE TO make_temperature_plot_terra: during scanner
        %%% change something was changed on putin for the temperature
        %%% configuration and the output file is different
        
        temp_trace = importdata(filename);
        
        if strcmp(scanner,'Terra_before_Sep20') %%% with old cable
            lines_to_use_idx = find(isnan(temp_trace.data(:,1)) | temp_trace.data(:,1) < 0); %%% the real trace is the one that only has one column in the stupid csv file, the part before is nonsense
            if lines_to_use_idx(end) ~= length(temp_trace.data) %%% workaround for one weird file
                lines_to_use_idx = lines_to_use_idx(end)+1 :  length(temp_trace.data);
            end
            temperature_trace = temp_trace.data(lines_to_use_idx,2); 
            temperature_trace_times = temp_trace.rowheaders(lines_to_use_idx);
        elseif strcmp(scanner,'Terra') %%% with new cable
            lines_to_use_idx = find(isnan(temp_trace.data(:,2)) | temp_trace.data(:,2) < 0);
            temperature_trace = temp_trace.data(lines_to_use_idx,1); 
            temperature_trace_times = temp_trace.textdata(lines_to_use_idx);
        elseif strcmp(scanner,'Weird')
            lines_to_use_idx = find(~isnan(temp_trace.data(:,1)) & temp_trace.data(:,1) > 0);
            temperature_trace = temp_trace.data(lines_to_use_idx,1); 
            temperature_trace_times = temp_trace.textdata(lines_to_use_idx);
        else
            lines_to_use_idx = find(isnan(temp_trace.data(:,2))); %%% the real trace is the one that only has one column in the stupid csv file, the part before is nonsense
            if lines_to_use_idx(end) ~= length(temp_trace.data) %%% workaround for one weird file
                lines_to_use_idx = lines_to_use_idx(end)+1 :  length(temp_trace.data);
            end
            temperature_trace = temp_trace.data(lines_to_use_idx,1);
            temperature_trace_times = temp_trace.rowheaders(lines_to_use_idx);
        end
        
        %%% convert temperature times to numbers
        midnight_check = 0;
        for t = 1:length(temperature_trace_times)
            curr_time = temperature_trace_times{t};
            if strcmp(curr_time,'12:00:00 AM')
                midnight_check = midnight_check + 1;
                display('ghost hour')
            end
            temptimenum(t) = datenum(curr_time) + midnight_check; %%% after midnight add 1 to daycount; to reverse: datetime(temptimenum(1),'ConvertFrom','datenum')
        end
        
        if strcmp(scanner, 'Magnetom')
            temptimenum = temptimenum + 1.5/24;  %%% account for putin, who was 1h30 behind back then
        end
        %%% convert seriesinfo time to numbers
        midnight_check = 0;
        for s = 1:length(seriesinfo)
            if length(seriesinfo(s).time) > 0
                try
                    seriestimenum(s) = datenum(seriesinfo(s).time) + midnight_check;
                catch
                    seriestimenum(s) = seriesinfo(s).time; %%% for case where file was again saved in another weird format
                end
                if s>1 && (seriestimenum(s-1)>seriestimenum(s) ) %%% if the previous sequence has higher timestamp, midnight passed
                    midnight_check = midnight_check + 1;
                    seriestimenum(s) = seriestimenum(s) + 1;
                end
                seriesname{s} = seriesinfo(s).name;
            else
               seriestimenum(s) = NaN;
               seriesname{s} = 'killme';
            end
        end
        seriestimenum = seriestimenum(~isnan(seriestimenum));
        seriesname = seriesname(~isnan(seriestimenum))
        
        %%% get temp for each sequence
        for s = 1:length(seriesname)
            temp_idx = find(temptimenum == seriestimenum(s));
            if length(temp_idx) == 0
               temp_idx = find(abs(temptimenum - seriestimenum(s))<10^-5); 
            end
            temp_per_seq(s).seq = seriesname{s};
            temp_per_seq(s).temp = temperature_trace(temp_idx);
        end
        %%% save in .mat file a vector with temperature for each sequence
        save(outputfilename_temp,'temp_per_seq');

        %%% find unique sequences
        idx = find([1,diff(seriestimenum)]>0);
        seriestimenum = seriestimenum(idx);
        seriesname = seriesname(idx);
        
        %%% make plot with 2 axes
        f = figure(1);
        %set(gcf,'Visible','off');
        set(groot,'defaultAxesTickLabelInterpreter','none');
        set(f,'position',[0 0 2000 800]);
        plot(temptimenum,temperature_trace);
        ylim([0 50])
        ylabel('Degrees Celsius');
        set(gca,'Position', get(gca,'Position').*[1 4 1 .5]) %%% make more space for labels
        ax1 = gca;
        ax1_pos = ax1.Position;
        set(ax1,'XTick',seriestimenum, 'XTickLabels',seriesname,'XTickLabelRotation',60);
        set(ax1,'FontSize',8);
        xlabel('Sequence');
        ax2 = axes('Position',ax1.Position,'XAxisLocation','top','YAxisLocation','right','Color','none', 'XLim',get(ax1,'XLim'),'YLim',get(ax1,'YLim'));
        set(ax2,'XTick',temptimenum(1:5000:length(lines_to_use_idx)), 'XTickLabels',temperature_trace_times(1:5000:length(lines_to_use_idx)),'XTickLabelRotation',60);
        set(ax2,'FontSize',8);
        xlabel('Actual time');
        ylabel('Degrees Celsius');

        saveas(f, outputfilename_fig)
        close all
end