clear all, close all, clc

Fs      = 250;
nbchans = 21;
nbpts   = 20;
win     = 5 * Fs;       % 5-sec windows                           
step    = 2.5 * Fs;     % 2.5-sec steps                         
nbwins  = (30*60*Fs-step) / step;             
TW      = 2.5;
FreqRes = 2 * TW / (win/Fs);   

% Post-COVID insomnia
load COVID_notch.mat

COVID1 = COVID_notch;     % time * chan * pt

t = 1:win;
COVID2 = NaN(win, nbwins, nbchans, nbpts);
for wini = 1:nbwins
    COVID2(:, wini, :, :) = COVID1(t, :, :);
    t = t + step;
end

% Multitaper power spectral density estimate
for chani = 1:nbchans
    for pi = 1:nbpts
        [pxx_COVID(:,:,chani,pi), f] = pmtm(COVID2(:,:,chani,pi)-mean(COVID2(:,:,chani,pi)), TW, win, Fs);
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

Pwr_COVID = cat(3, f1, f2, f3, f4, f5);



% Primary insomnia
load Insom_notch.mat

Insom1 = Insom_notch;

t = 1:win;
Insom2 = NaN(win, nbwins, nbchans, nbpts);
for wini = 1:nbwins
    Insom2(:, wini, :, :) = Insom1(t, :, :);
    t = t + step;
end

for chani = 1:nbchans
    for pi = 1:nbpts
        [pxx_Insom(:,:,chani,pi), f] = pmtm(Insom2(:,:,chani,pi)-mean(Insom2(:,:,chani,pi)), TW, win, Fs);
    end
end

f1 = squeeze( mean(mean(pxx_Insom(fidx1(1):fidx1(2),:,:,:), 1), 2) );
f2 = squeeze( mean(mean(pxx_Insom(fidx2(1):fidx2(2),:,:,:), 1), 2) );
f3 = squeeze( mean(mean(pxx_Insom(fidx3(1):fidx3(2),:,:,:), 1), 2) );
f4 = squeeze( mean(mean(pxx_Insom(fidx4(1):fidx4(2),:,:,:), 1), 2) );
f5 = squeeze( mean(mean(pxx_Insom(fidx5(1):fidx5(2),:,:,:), 1), 2) );

Pwr_Insom = cat(3, f1, f2, f3, f4, f5);



% Healthy controls
load NC_notch.mat

NC1 = NC_notch;

t = 1:win;
Insom2 = NaN(win, nbwins, nbchans, nbpts);
for wini = 1:nbwins
    NC2(:, wini, :, :) = NC1(t, :, :);
    t = t + step;
end

for chani = 1:nbchans
    for pi = 1:nbpts
        [pxx_NC(:,:,chani,pi), f] = pmtm(NC2(:,:,chani,pi)-mean(NC2(:,:,chani,pi)), TW, win, Fs);
    end
end

f1 = squeeze( mean(mean(pxx_NC(fidx1(1):fidx1(2),:,:,:), 1), 2) );
f2 = squeeze( mean(mean(pxx_NC(fidx2(1):fidx2(2),:,:,:), 1), 2) );
f3 = squeeze( mean(mean(pxx_NC(fidx3(1):fidx3(2),:,:,:), 1), 2) );
f4 = squeeze( mean(mean(pxx_NC(fidx4(1):fidx4(2),:,:,:), 1), 2) );
f5 = squeeze( mean(mean(pxx_NC(fidx5(1):fidx5(2),:,:,:), 1), 2) );

Pwr_NC = cat(3, f1, f2, f3, f4, f5);