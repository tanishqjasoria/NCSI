% model fiber parameters
clear all; 

CF_1    = 500; % CF in Hz; 
CF_2    = 4000; % CF in Hz;

stimdb_start = -10;
stimdb_stop  =  80;
stimdb_step  =  10;

freq_start   =  62.5;
freq_count   =  9;
freq_step    =  1.0/8;

cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
fiberType = 3; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse

% stimulus parameters
F0 = CF_1;     % stimulus frequency in Hz
Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = 200e-3;  % stimulus duration in seconds
rt = 10e-3;   % rise/fall time in seconds
stimdb = 10; % stimulus intensity in dB SPL

% PSTH parameters
nrep = 50               % number of stimulus repetitions (e.g., 50);
psthbinwidth = 0.1e-3; % binwidth in seconds;

% Experimental Parameters
stimdb_range    = stimdb_start:stimdb_step:stimdb_stop;
frequency_range = freq_start*2.^(0:freq_step:freq_count);


data_1 = zeros(length(frequency_range), length(stimdb_range));

for frequency_index=1:length(frequency_range)
    data_intensity = zeros(1, length(stimdb_range));
    frequency = frequency_range(frequency_index);
    for intensity_index=1:length(stimdb_range)
        intensity = stimdb_range(intensity_index);
        
        pin = get_stim(frequency, Fs, T, rt, intensity);
        [synout, psth] = ANModel(nrep, pin, CF_1, Fs, T, cohc, cihc, fiberType, implnt, psthbinwidth);
        
        data_intensity(1, intensity_index) = sum(psth);
    end
    data_1(frequency_index, :) = data_intensity;
end

data_2 = zeros(length(frequency_range), length(stimdb_range));

for frequency_index=1:length(frequency_range)
    data_intensity = zeros(1, length(stimdb_range));
    frequency = frequency_range(frequency_index);
    for intensity_index=1:length(stimdb_range)
        intensity = stimdb_range(intensity_index);
        
        pin = get_stim(frequency, Fs, T, rt, intensity);
        [synout, psth] = ANModel(nrep, pin, CF_2, Fs, T, cohc, cihc, fiberType, implnt, psthbinwidth);
        
        data_intensity(1, intensity_index) = sum(psth);
    end
    data_2(frequency_index, :) = data_intensity;
end
figure
hold on;

xtick = 62.5*2.^(0:(1.0/8):9);
set(gca,'XTick',xtick);
contour(gca, data_1');
 
figure
hold on;

xtick = 62.5*2.^(0:(1.0/8):9);
set(gca,'XTick',xtick)
%imagesc(gca, [frequency_range(1), frequency_range(length(frequency_range))], [stimdb_range(1) , stimdb_range(length(stimdb_range))] , data_1');
contour(gca, data_2')
%imagesc('XData',freqRange,'YData',intensityRange,'CData',experimentData')
 
% figure;
% for i=1:length(data_1)
%     hold on;
%     plot(stimdb_range, data_1(i, :));
% end

data_intensity_CF1 = zeros(1, length(stimdb_range));
for intensity_index=1:length(stimdb_range)
    intensity = stimdb_range(intensity_index);

    pin = get_stim(CF_1, Fs, T, rt, intensity);
    [synout, psth] = ANModel(nrep, pin, CF_1, Fs, T, cohc, cihc, fiberType, implnt, psthbinwidth);

    data_intensity_CF1(1, intensity_index) = sum(psth);
end

figure;
hold on;
plot(stimdb_range, data_intensity_CF1,'DisplayName','500');

data_intensity_CF2 = zeros(1, length(stimdb_range));
for intensity_index=1:length(stimdb_range)
    intensity = stimdb_range(intensity_index);

    pin = get_stim(CF_2, Fs, T, rt, intensity);
    [synout, psth] = ANModel(nrep, pin, CF_2, Fs, T, cohc, cihc, fiberType, implnt, psthbinwidth);

    data_intensity_CF2(1, intensity_index) = sum(psth);
end

plot(stimdb_range, data_intensity_CF2,'DisplayName','4000');
legend()

function [psthtime, Psth] = ANModel(nrep, pin, CF, Fs , T, cohc, cihc, fibreType, implnt, psthbinwidth)
    vihc = catmodel_IHC(pin,CF,nrep,1.0/Fs,T*2,cohc,cihc); 
    [synout,psth] = catmodel_Synapse(vihc,CF,nrep,1.0/Fs,fibreType, implnt);
    disp(length(psth))
    timeout = (1:length(psth))*1/Fs;
    psthbins = round(psthbinwidth*Fs);  % number of psth bins per psth bin
    psthtime = timeout(1:psthbins:end); % time vector for psth
    pr = sum(reshape(psth,psthbins,length(psth)/psthbins))/nrep; % pr of spike in each bin
    Psth = pr/psthbinwidth; % psth in units of spikes/s
end


