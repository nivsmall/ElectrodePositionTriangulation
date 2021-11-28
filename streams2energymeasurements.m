function [M] = streams2energymeasurements(parentDirSync, parentDirStream, len, harmonicCycleTime, sampleRate)
% fetching sine wave acttuation and measurements from continuous stream & sync
% performing integral over sine cycle to get total energy measurement
%   data of each trial held at different file
%   harmonicCycleTime - [seconds]
%   sampleRate - [Hz]
    streamFiles = dir(parentDirStream);
    syncFiles = dir(parentDirSync);
    n = length(streamFiles);
    energy = nan(len, len, n-2);
    cnt=1;
    for f = 1:n
        if streamFiles(f).isdir
            continue
        end
        syncfile = strcat(parentDirSync, syncFiles(f).name);
        streamfile = strcat(parentDirStream, streamFiles(f).name);
        mask = load(syncfile).mask;
        stream = load(streamfile).M;
        data = stream.*mask;
        
        i=1;
        while mask(i) == 0
            i=i+1;
        end
        firstSync = i;
        while mask(i) == 1
            i=i+1;
        end
        samplesInCycle = i-firstSync;
        while mask(i) == 0
            i=i+1;
        end
        secondSync = i;
        fullCycleTime = secondSync - firstSync;
        
        %energy per cycle:
        %samplesInCycle = harmonicCycleTime*sampleRate;
        cycles = ceil(sum(mask)/samplesInCycle); % #cycles = number of exciting electrodes
        
        for i = 0:cycles-1
            if i >= len
                break;
            end
            fullCycleMeasurment = data(:, firstSync+i*fullCycleTime : firstSync + i*fullCycleTime + samplesInCycle);
            energy(:, i+1, cnt) = sum(abs(fullCycleMeasurment), 2);
        end
        
        cnt=cnt+1;
    end
        
M = mean(energy, 3);

end

