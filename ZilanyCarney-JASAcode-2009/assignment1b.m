

% model fiber parameters
clear all; 
CF    = 500; % CF in Hz;   
cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
fiberType = 3; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse

% stimulus parameters
% characterstics frequency range must be within 80Hz - 40kHz (model
% restrictions)
frequency_range = 80*2.^(0:1/8:7); 
Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
rt = 10e-3;   % rise/fall time in seconds
stimdb_range = -20:10:80;

% PSTH parameters
nrep = 1               % number of stimulus repetitions (e.g., 50);
psthbinwidth = 0.1e-3; % binwidth in seconds;

[audio, Fs_audio] = audioread('fivewo.wav');
audioinfo('fivewo.wav')

T = length(audio)*(1/Fs_audio);

ah_audio = audio(110001:1:120000)';

rms_ah = rms(ah_audio);
dB_spl = 20 * log10(rms_ah/20*10^(-6));

mxpts = length(ah_audio);
irpts = double(int64(rt*Fs_audio));


ah_audio(1:irpts) = ah_audio(1:irpts).*(0:(irpts-1))/irpts; 
ah_audio((mxpts-irpts):mxpts) = ah_audio((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;

data = ones(1, length(stimdb_range));
%%
figure;
hold on;

for j=1:length(stimdb_range)
    intensity = stimdb_range(j);
    pin = ah_audio * sqrt(2)*20e-6*10^(intensity/20.0)/rms_ah;
    [synout, psth] = ANModel(nrep, pin, CF, Fs, length(pin)/Fs, cohc, cihc, fiberType,implnt, psthbinwidth); 
    data(1,j) = sum(psth);
end


plot(stimdb_range, data,'DisplayName','ah')

for j=1:length(stimdb_range)
    intensity = stimdb_range(j);
    pin = get_stim(CF, Fs, length(ah_audio)/Fs, rt, intensity);
    [synout, psth] = ANModel(nrep, pin, CF, Fs, length(pin)/Fs, cohc, cihc, fiberType,implnt, psthbinwidth); 
    data(1,j) = sum(psth);
end
plot(stimdb_range, data,'DisplayName','BF')
legend()

window_size = 25.6e-3 * Fs;

overlap_percent = 50/100;
overlap_size = window_size * overlap_percent;

figure;
spectrogram(audio, window_size, overlap_size, window_size, Fs, 'yaxis');
spec_axis = gca;
spec_axis.YScale = 'log';


F_Th = 5;
F_DR = 35;
F_SL = 75;

sound_th = sqrt(2)*20e-6*10^(F_Th/20.0)* audio / rms_ah;
sound_dr = sqrt(2)*20e-6*10^(F_DR/20.0)* audio / rms_ah;
sound_sl = sqrt(2)*20e-6*10^(F_SL/20.0)* audio / rms_ah;

%%
nrep = 50
figure;

windows_time = [4e-3, 8e-3, 16e-3, 32e-3, 64e-3, 128e-3];
overlap_percent = 50/100;
for time_index=1:length(windows_time)
    
    window_size = windows_time(time_index) * Fs;
    overlap_size = window_size * overlap_percent;
    
    windows_size = window_size - overlap_size;
    
    [psthtime, psth] = ANModel(nrep, sound_th', frequency_range(1), Fs, length(sound_th)/Fs, cohc, cihc, fiberType,implnt, psthbinwidth);
    length_psth_sum = length(1:windows_size:length(psth));
    
    data = zeros(length(frequency_range), length_psth_sum);
    
    for frequency_index=1:length(frequency_range)
        CF = frequency_range(frequency_index);
        [psthtime, psth] = ANModel(nrep, sound_th', CF, Fs, length(sound_th)/Fs, cohc, cihc, fiberType,implnt, psthbinwidth);
        psth_windows = window_psth(overlap_size, window_size, psth);
        data( frequency_index, :) = psth_windows;
        display(length(psth_windows));
    end
    
    subplot(3,2, time_index);
    title( strcat(int2str(window_size/100),'ms'))
    ytick = 80*2.^(0:1/8:7);
    set(gca,'YTick',ytick);
    imagesc(gca, 'XData', [1, length(psth_windows)], 'YData', [frequency_range(1), frequency_range(length(frequency_range))], 'CData', data);

end

%%
figure;

windows_time = [4e-3, 8e-3, 16e-3, 32e-3, 64e-3, 128e-3];
overlap_percent = 50/100;
for time_index=1:length(windows_time)
    
    window_size = windows_time(time_index) * Fs;
    overlap_size = window_size * overlap_percent;
    
    windows_size = window_size - overlap_size;
    
    [psthtime, psth] = ANModel(nrep, sound_dr', frequency_range(1), Fs, length(sound_dr)/Fs, cohc, cihc, fiberType,implnt, psthbinwidth);
    length_psth_sum = length(1:windows_size:length(psth));
    
    data = zeros(length(frequency_range), length_psth_sum);
    
    for frequency_index=1:length(frequency_range)
        CF = frequency_range(frequency_index);
        [psthtime, psth] = ANModel(nrep, sound_dr', CF, Fs, length(sound_dr)/Fs, cohc, cihc, fiberType,implnt, psthbinwidth);
        psth_windows = window_psth(overlap_size, window_size, psth);
        data( frequency_index, :) = psth_windows;
        display(length(psth_windows));
    end
    
    subplot(3,2, time_index);
    title( strcat(int2str(window_size/100),'ms'))
    ytick = 80*2.^(0:1/8:7);
    set(gca,'YTick',ytick);
    imagesc(gca, 'XData', [1, length(psth_windows)], 'YData', [frequency_range(1), frequency_range(length(frequency_range))], 'CData', data);

end

%%
figure;

windows_time = [4e-3, 8e-3, 16e-3, 32e-3, 64e-3, 128e-3];
overlap_percent = 50/100;
for time_index=1:length(windows_time)
    
    window_size = windows_time(time_index) * Fs;
    overlap_size = window_size * overlap_percent;
    
    windows_size = window_size - overlap_size;
    
    [psth, psth] = ANModel(nrep, sound_sl', frequency_range(1), Fs, length(sound_sl)/Fs, cohc, cihc, fiberType,implnt, psthbinwidth);
    length_psth_sum = length(1:windows_size:length(psth));
    
    data = zeros(length(frequency_range), length_psth_sum);
    
    for frequency_index=1:length(frequency_range)
        CF = frequency_range(frequency_index);
        [time, psth] = ANModel(nrep, sound_sl', CF, Fs, length(sound_sl)/Fs, cohc, cihc, fiberType,implnt, psthbinwidth);
        psth_windows = window_psth(overlap_size, window_size, psth);
        data( frequency_index, :) = psth_windows;
        display(length(psth_windows));
    end
    
    subplot(3,2, time_index);
    title( strcat(int2str(window_size/100),'ms'))
    ytick = 80*2.^(0:1/8:7);
    set(gca,'YTick',ytick);
    imagesc(gca, 'XData', [1, length(psth_windows)], 'YData', [frequency_range(1), frequency_range(length(frequency_range))], 'CData', data);

end

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

function window_psth_sum = window_psth(overlap_size, binwidth, psth)
%WINDOW_PSTH Summary of this function goes here
%   Detailed explanation goes here
    windows_size = binwidth - overlap_size;
    window_starts = 1:windows_size:length(psth);
    
    window_psth_sum = zeros(1, length(window_starts));
    
    for i=1:length(window_starts)
        window_start = window_starts(i);
        window_end = window_start + binwidth - 1;
        
        if window_end <= length(psth)
            window_psth_sum(i) = sum(psth(window_start: window_end));
        else
            window_psth_sum(i) = sum(psth(window_start:end));
        end
    end
end