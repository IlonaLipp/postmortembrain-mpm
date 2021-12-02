function ReadCoilTempFile(logfilename, dateofinterest, figurefilename)

    %%% this was created by Evgeniya to plot gradient coil temperatures, based
    %%% on the log file saved from siemens. it reads information from all
    %%% sensors and plots them as a function of time

    %% read CoilTempFile
    fid = fopen(logfilename);
    tline = fgetl(fid);
    i=1;
    while ischar(tline)
        %disp(tline)
        tline = fgetl(fid);
       if (strfind(tline,dateofinterest)>0)
            if (strfind(tline,'GC1:')>0)
            time=tline((strfind(tline, dateofinterest)+11):(strfind(tline, dateofinterest)+25));
            [Y, M, D, H, MN, S] = datevec(time);
            t(i)=H*60+MN+S/60;
            temp(1,i)=str2num(tline((strfind(tline,'GC1:')+5):(strfind(tline,'GC1:')+9)));
            temp(2,i)=str2num(tline((strfind(tline,'GC2:')+5):(strfind(tline,'GC2:')+9)));
            temp(3,i)=str2num(tline((strfind(tline,'GC3:')+5):(strfind(tline,'GC3:')+9)));
            temp(4,i)=str2num(tline((strfind(tline,'GC4:')+5):(strfind(tline,'GC4:')+9)));
            temp(5,i)=str2num(tline((strfind(tline,'GC5:')+5):(strfind(tline,'GC5:')+9)));
            temp(6,i)=str2num(tline((strfind(tline,'GC6:')+5):(strfind(tline,'GC6:')+9)));
            temp(7,i)=str2num(tline((strfind(tline,'GC7:')+5):(strfind(tline,'GC7:')+9)));
            temp(8,i)=str2num(tline((strfind(tline,'GC8:')+5):(strfind(tline,'GC8:')+9)));
            temp(9,i)=str2num(tline((strfind(tline,'GC9:')+5):(strfind(tline,'GC9:')+9)));
            temp(10,i)=str2num(tline((strfind(tline,'GC10:')+5):(strfind(tline,'GC10:')+9)));
            temp(11,i)=str2num(tline((strfind(tline,'GC11:')+5):(strfind(tline,'GC11:')+9)));
            temp(12,i)=str2num(tline((strfind(tline,'GC12:')+5):(strfind(tline,'GC12:')+9)));
            i=i+1;
           end;
       end;
    end
    fclose(fid);
    
    f=figure()
    a=plot(t,temp(1,:),'-or', t, temp(2,:),'-og', t, temp(3,:),'-ob', t, temp(4,:),'-oc',t, temp(5,:),'-om', t, temp(6,:), '-oy',...
         t,temp(7,:),'-ok', t, temp(8,:),'-*r', t, temp(9,:),'-*g', t, temp(10,:),'-*b',t, temp(11,:),'-*c', t, temp(12,:), '-*m', t, 50*ones(length(t)),'--r');
     legend(a,{'GC1', 'GC2', 'GC3', 'GC4', 'GC5', 'GC6', 'GC7', 'GC8', 'GC9', 'GC10', 'GC11', 'GC12'},'Location','Northwest');
     axis([min(t) max(t) 15 70]);

     title(strcat('Gradient Coil Temperature at 7T Terra on:', '    ', dateofinterest));
     xlabel('Time, min') ;
     ylabel('Temp, �C') 
     saveas(f,figurefilename)
end
