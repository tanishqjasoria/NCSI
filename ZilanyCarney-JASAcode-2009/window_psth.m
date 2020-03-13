function window_psth_sum = window_psth(overlap_samples, binwidth_samples, psth)
%WINDOW_PSTH Summary of this function goes here
%   Detailed explanation goes here
    windows_size = binwidth_samples - overlap_samples;
    window_starts = 1:windows_size:length(psth);
    
    window_psth_sum = zeros(1, length(window_starts));
    
    for i=1:length(window_starts)
        window_start = window_starts(i);
        window_end = window_start + binwidth_samples - 1;
        
        if window_end <= length(psth)
            window_psth_sum(i) = sum(psth(window_start: window_end));
        else
            window_psth_sum(i) = sum(psth(window_start:end));
        end
    end
end
    

