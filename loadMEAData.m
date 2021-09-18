
read_data = false;
pulseData = false;
num_trials = 10;

%% Read the data:
%%%% Data Flow prior MATLAB: MCRack -> .cmd file -> MCDataManager -> .h5 file %%%%
samplerate = 1/10000;       % 10[kHz]
if read_data
    if pulseData
        load('continuousStream51to60.mat');
        load('continuousStream1to50.mat');
        load('triggers_51to60.mat');
        load('triggers_1to50.mat');
        load('ChannelDataTimeStamps_51to60.mat');
        load('ChannelDataTimeStamps_1to50.mat');
    end
else
    if pulseData
        % configure the segment of data to read:
        cfg_1to50=[];
        cfg_1to50.window=[55496  15008090].*samplerate; % data slicing w.r.t time [s] [5549.7, 1999894.1] / 2031874 / 5026385 / 10037200]
        cfg_51to60=[];
        cfg_51to60.window=[42200 3016750].*samplerate;  % 1st trigger @ [4220.1, 401491.9] 421506 / 102040.8 / 2218212]
        data_path_1to50 = 'C:/Users/nivsm/Documents/FinalProject/MEA_experiments/MultichannelDataManager/50_first_excites_approx20trials0001.h5';
        data_path_51to60 = 'C:/Users/nivsm/Documents/FinalProject/MEA_experiments/MultichannelDataManager/51to60_excites_approx20trials.h5';
        data_1to50 = McsHDF5.McsData(data_path_1to50);
        data_51to60 = McsHDF5.McsData(data_path_51to60);
        % triggers:                 1 X #excitations
        triggers_51to60 = data_51to60.Recording{1, 1}.EventStream{1, 1}.readPartialEventData(cfg_51to60).Events{1, 1};
        %save('triggers_51to60.mat', 'triggers_51to60');
        triggers_1to50 = data_1to50.Recording{1, 1}.EventStream{1, 1}.readPartialEventData(cfg_1to50).Events{1, 1};
        %save('triggers_1to50.mat', 'triggers_1to50');
        % continuous data stream:   #channels X #samples
        M_51to60 = data_51to60.Recording{1, 1}.AnalogStream{1, 2}.readPartialChannelData(cfg_51to60);
        M_51to60Cont =M_51to60.ChannelData;
        ChannelDataTimeStamps_51to60 = M_51to60.ChannelDataTimeStamps;
        save('ChannelDataTimeStamps_51to60.mat', 'ChannelDataTimeStamps_51to60');
        save('continuousStream51to60.mat', 'M_51to60Cont')
        M_1to50 = data_1to50.Recording{1, 1}.AnalogStream{1, 2}.readPartialChannelData(cfg_1to50);
        M_1to50Cont =M_1to50.ChannelData;
        ChannelDataTimeStamps_1to50 = M_1to50.ChannelDataTimeStamps;
        save('ChannelDataTimeStamps_1to50.mat', 'ChannelDataTimeStamps_1to50');
        save('continuousStream1to50.mat', 'M_1to50Cont');
    else
        
        parentDir = 'C:/Users/nivsm/Documents/FinalProject/MEA_experiments/MultichannelDataManager/harmonicStim/';
        dataFiles = dir(parentDir);
        n = length(dataFiles);
        cnt = 1;
        for f = 1:n
            if dataFiles(f).isdir
                continue
            end
            filePath = strcat(parentDir, dataFiles(f).name);
            data = McsHDF5.McsData(filePath);
            % triggers:                 1 X #excitations
            %triggers = data.Recording{1, 1}.EventStream{1, 1}.readPartialEventData(CFG_HERE).Events{1, 1};
            % digital sync (for timing stimulation):
            mask = data.Recording{1, 1}.AnalogStream{1, 1}.ChannelData + 32768;
            % continuous data stream:   #channels X #samples
            M = data.Recording{1, 1}.AnalogStream{1, 2}.ChannelData;
            %M_Cont = M.ChannelData;
            %ChannelDataTimeStamps = M.ChannelDataTimeStamps;
            nameToWrite = strcat('syncData/stimulusSyncMask', string(cnt), '.mat');
            save(nameToWrite, 'mask'); %sync with id
            nameToWrite = strcat('measData/dataStream', string(cnt), '.mat');
            save(nameToWrite, 'M'); %stream with id
            
            cnt=cnt+1;
        end
    
    end
end

if pulseData
    % extract pulses (peaks), rearrange data:
    trig_indices_1to50 = trigger2triggerIndices(triggers_1to50, ChannelDataTimeStamps_1to50);
    trig_indices_51to60 = trigger2triggerIndices(triggers_51to60, ChannelDataTimeStamps_51to60);
    M_1to50 = continuousStream2triggered(M_1to50Cont, trig_indices_1to50, 50, num_trials);
    M_51to60 = continuousStream2triggered(M_51to60Cont, trig_indices_51to60, 10, num_trials);
else
    % extract recorded energy (over time), rearrange data:
    Pot_Mat = streams2energymeasurements('syncData/', 'measData/', 60, 1, 10*10^3);
    save('Pot_Mat_harmonicMeas.mat', 'Pot_Mat');
    return;
end

% Potential Matrix - coloumn # is the exciting electrode, row # is the
% measuring electrode
Pot_Mat = abs(horzcat(M_1to50, M_51to60));
save('Pot_Mat_experiment_24_2_2021_X2.mat', 'Pot_Mat')

return


