function [synout,psth] = ANModel(nrep, pin, CF, Fs , T, cohc, cihc, fibreType, implnt)
    vihc = catmodel_IHC(pin,CF,nrep,1.0/Fs,T*2,cohc,cihc); 
    [synout,psth] = catmodel_Synapse(vihc,CF,nrep,1.0/Fs,fibreType, implnt); 
end