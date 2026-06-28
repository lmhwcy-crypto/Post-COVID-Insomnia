clear all, close all, clc

% Post-COVID insomnia
% Load your EEG matrix, sleep or waking EEG: time points * channels * patients 
load COVID_notch.mat

Fs      = 250;          % Your sampling rate
nbchans = size(COVID_notch, 2);
nbpts   = size(COVID_notch, 3);
win     = 5 * Fs;       % 5-sec windows                           
step    = 2.5 * Fs;     % 2.5-sec steps                         
nbwins  = (30*60*Fs-step) / step;             
TW      = 2.5;
FreqRes = 2 * TW / (win/Fs);   

t = 1:win;
COVID = NaN(win, nbwins, nbchans, nbpts);
for wini = 1:nbwins
    COVID(:, wini, :, :) = COVID_notch(t, :, :);
    t = t + step;
end

% Multitaper power spectral density estimate
for chani = 1:nbchans
    for pi = 1:nbpts
        [pxx_COVID(:,:,chani,pi), f] = pmtm(COVID(:,:,chani,pi)-mean(COVID(:,:,chani,pi)), TW, win, Fs);
        % Frequency * time windows * channels * patients
    end
end

% Find frequencies for delta, theta, alpha, beta, and gamma
fidx1 = dsearchn(f, [0.6; 3.8]);       
fidx2 = dsearchn(f, [4; 7.8]);   
fidx3 = dsearchn(f, [8; 13]);   
fidx4 = dsearchn(f, [13.2; 30]);   
fidx5 = dsearchn(f, [30.2; 80]);   

f1 = squeeze( mean(mean(pxx_COVID(fidx1(1):fidx1(2),:,:,:), 1), 2) );
f2 = squeeze( mean(mean(pxx_COVID(fidx2(1):fidx2(2),:,:,:), 1), 2) );
f3 = squeeze( mean(mean(pxx_COVID(fidx3(1):fidx3(2),:,:,:), 1), 2) );
f4 = squeeze( mean(mean(pxx_COVID(fidx4(1):fidx4(2),:,:,:), 1), 2) );
f5 = squeeze( mean(mean(pxx_COVID(fidx5(1):fidx5(2),:,:,:), 1), 2) );

% Get the power density matrix: Channels * patients * frequencies
Pwr_COVID = cat(3, f1, f2, f3, f4, f5);



% Primary insomnia
% Load your EEG matrix: time points * channels * patients 
load Insom_notch.mat

t = 1:win;
Insom = NaN(win, nbwins, nbchans, nbpts);
for wini = 1:nbwins
    Insom(:, wini, :, :) = Insom_notch(t, :, :);
    t = t + step;
end

% Multitaper power spectral density estimate
for chani = 1:nbchans
    for pi = 1:nbpts
        [pxx_Insom(:,:,chani,pi), f] = pmtm(Insom(:,:,chani,pi)-mean(Insom(:,:,chani,pi)), TW, win, Fs);
        % Frequency * time windows * channels * patients
    end
end

f1 = squeeze( mean(mean(pxx_Insom(fidx1(1):fidx1(2),:,:,:), 1), 2) );
f2 = squeeze( mean(mean(pxx_Insom(fidx2(1):fidx2(2),:,:,:), 1), 2) );
f3 = squeeze( mean(mean(pxx_Insom(fidx3(1):fidx3(2),:,:,:), 1), 2) );
f4 = squeeze( mean(mean(pxx_Insom(fidx4(1):fidx4(2),:,:,:), 1), 2) );
f5 = squeeze( mean(mean(pxx_Insom(fidx5(1):fidx5(2),:,:,:), 1), 2) );

% Get the power density matrix: Channels * patients * frequencies
Pwr_Insom = cat(3, f1, f2, f3, f4, f5);



% Healthy controls
% Load your EEG matrix: time points * channels * patients 
load NC_notch.mat

t = 1:win;
NC = NaN(win, nbwins, nbchans, nbpts);
for wini = 1:nbwins
    NC(:, wini, :, :) = NC_notch(t, :, :);
    t = t + step;
end

% Multitaper power spectral density estimate
for chani = 1:nbchans
    for pi = 1:nbpts
        [pxx_NC(:,:,chani,pi), f] = pmtm(NC(:,:,chani,pi)-mean(NC(:,:,chani,pi)), TW, win, Fs);
        % Frequency * time windows * channels * patients
    end
end

f1 = squeeze( mean(mean(pxx_NC(fidx1(1):fidx1(2),:,:,:), 1), 2) );
f2 = squeeze( mean(mean(pxx_NC(fidx2(1):fidx2(2),:,:,:), 1), 2) );
f3 = squeeze( mean(mean(pxx_NC(fidx3(1):fidx3(2),:,:,:), 1), 2) );
f4 = squeeze( mean(mean(pxx_NC(fidx4(1):fidx4(2),:,:,:), 1), 2) );
f5 = squeeze( mean(mean(pxx_NC(fidx5(1):fidx5(2),:,:,:), 1), 2) );

% Get the power density matrix: Channels * patients * frequencies
Pwr_NC = cat(3, f1, f2, f3, f4, f5);
