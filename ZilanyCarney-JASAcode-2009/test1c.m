%
clear all;
CF = 500; % CF in Hz;
cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
fiberType = 2; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
%

%
% stimulus parameters
frequency_range = 80*2.^(0:1/8:7);
frequency_range_mark = 125*2.^(0:1/2:4);
stimdb_range = -50:10:200;
rt = 10e-3;   % rise/fall time in seconds

% PSTH parameters
nrep = 1;               % number of stimulus repetitions (e.g., 50);
psthbinwidth = 0.01e-3; % binwidth in seconds;

%%

[audio, Fs_audio] = audioread('fivewo.wav');
audioinfo('fivewo.wav')

Fs = 500e3;
PSTH = zeros(50, 62500);
for i=1:50
    [psthtime, psth] = ANModel(nrep, audio', frequency_range(1), Fs, length(audio)/Fs, cohc, cihc, fiberType, implnt, psthbinwidth);
    PSTH(i,:) = psth;
end

%%
% Fs = 100e3;
% window_samples = 0.01e-3*Fs;
% overlap_samples = floor(window_samples/2);
% 
% [synout, psth] = ANModel(nrep, audio', frequency_range(1), Fs, length(audio)/Fs, cohc, cihc, fiberType,implnt);     
% psth = window_psth(overlap_samples, window_samples, psth);
% data = zeros(length(frequency_range), length(psth));
% data_1 = zeros(length(frequency_range_1), length(psth));
% data(1,:) = psth;
% 
% %%
% % ANF model 1
% parfor i=2:length(frequency_range)
%     CF = frequency_range(i);
%     [synout, psth] = ANModel(nrep, audio', CF, Fs, length(audio)/Fs, cohc, cihc, fiberType,implnt);     
%     psth = window_psth(overlap_samples, window_samples, psth);
%     data(i,:) = psth;
% end
% %%
% figure;
% 
% imagesc('XData', [1,length(psth)], 'YData', [frequency_range(1), frequency_range(length(frequency_range))], 'CData', flipud(data))
% ytick = 80*2.^(0:1/8:7);
% set(gca,'YTick',ytick);
% hold on;
% %%
% % ANF model 1
% parfor i=1:length(frequency_range_1)
%     CF = frequency_range_1(i);
%     [synout, psth] = ANModel(nrep, audio', CF, Fs, length(audio)/Fs, cohc, cihc, fiberType,implnt);     
%     psth = window_psth(overlap_samples, window_samples, psth);
%     data_1(i,:) = psth;
% end

%% 


fftWindow = 12.8e-3*Fs;
fftoverlap = fftWindow/2;
indices = 1:(fftWindow-fftoverlap):length(data_1(1,:));

L = fftWindow;                
T = 1/Fs;             % Sampling period       
t = (0:L-1)*T;

L1 = length(data_1(1,:));

for i=1:length(frequency_range_mark)
    col = rand(1,3);
    for j=1:length(indices)
        if indices(j)+fftWindow-1 > L1
            Y = fft(data_1(i,indices(j):end));   
            L = length(data_1(i,indices(j):end));                
            T = 1/Fs;             % Sampling period       
            t = (0:L-1)*T;
        else
            Y = fft(data_1(i,indices(j):indices(j)+fftWindow-1));   
        end
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L;
        [argval, argmax] = max(P1);
        if f(argmax) < 20000 && f(argmax) ~= 0
            if j > 1
                plot((indices(j-1) + indices(j))/2, f(argmax),'*','color',col)
            else
                plot((indices(j))/2, f(argmax),'*','color',col)
            end
        end
        
    end
end

%% 