
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
nrep = 1;               % number of stimulus repetitions (e.g., 50);
psthbinwidth = 0.1e-3; % binwidth in seconds;

[audio, Fs] = audioread('fivewo.wav');
audioinfo('fivewo.wav')

T = length(audio)*(1/Fs);

ah_audio = audio(110000:1:120000)';

rms_ah = rms(ah_audio);
dB_spl = 20 * log10(rms_ah/20*10^(-6));

mxpts = length(ah_audio);
irpts = double(int64(rt*Fs));


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

%%

F_Th = 40;
F_DR = 100;
F_SL = 160;

sound_th = sqrt(2)*20e-6*10^(F_Th/20.0)* audio / rms_ah;
sound_dr = sqrt(2)*20e-6*10^(F_DR/20.0)* audio / rms_ah;
sound_sl = sqrt(2)*20e-6*10^(F_SL/20.0)* audio / rms_ah;

window_samples = 25.6e-3 * Fs;

overlap_percent = 50/100;
overlap_samples = window_samples * overlap_percent;

figure;
spectrogram(audio, window_samples, overlap_samples, window_samples, Fs);
ax = gca;
ax.XScale = 'log';



function [Psth, psthtime] = ANModel(nrep, pin, CF, Fs , T, cohc, cihc, fibreType, implnt, psthbinwidth)
    vihc = catmodel_IHC(pin,CF,nrep,1.0/Fs,T*2,cohc,cihc); 
    [synout,psth] = catmodel_Synapse(vihc,CF,nrep,1.0/Fs,fibreType, implnt);
    
    timeout = (1:length(psth))*1/Fs;
    psthbins = round(psthbinwidth*Fs);  % number of psth bins per psth bin
    psthtime = timeout(1:psthbins:end); % time vector for psth
    pr = sum(reshape(psth,psthbins,length(psth)/psthbins))/nrep; % pr of spike in each bin
    Psth = pr/psthbinwidth; % psth in units of spikes/s
end
