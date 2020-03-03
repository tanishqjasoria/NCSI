function pin = generateStimulus(F0,Fs, T, rt, stimdb)
    
    t = 0:1/Fs:T-1/Fs; % time vector
    mxpts = length(t);
    irpts = rt*Fs;

    pin = sqrt(2)*20e-6*10^(stimdb/20.0)*sin(2*pi*F0*t); % unramped stimulus
    pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts; 
    pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;

end

