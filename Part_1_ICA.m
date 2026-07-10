clear all, close all, clc

% Load your EEG matrix, sleep or waking EEG: record1 = channels * time points

% Filter EEG with a Notch filter at 60 Hz

srate = 250;
fcomb = [55 59 61 64];
mags = [1 0 1];
dev = [0.5 0.1 0.5];
[n,Wn,beta,ftype] = kaiserord(fcomb,mags,dev,srate);
hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');

record2 = FiltFiltM(hh, 1, record1');
record  = record2';

%%
%ICA

clear all, close all, clc

load ChanLabels_CGMH.mat     % Channel label vetor

% Load your EEG matrix, sleep or waking EEG: record = channels * time points
nbchans = 21;
data    = record(1:nbchans, :);   

weight = NaN(nbchans, nbchans);
sphere = NaN(nbchans, nbchans);
W      = NaN(nbchans, nbchans);
icaEEG = NaN(nbchans, length(data)); 
pvaf   = NaN(1, nbchans);   

[weight, sphere] = runica(data, 'verbose', 'off');
W = weight * sphere;
icaEEG = W * data;
winv = inv(W);
squaredata = sum(sum(data.^2));

for compi = 1:nbchans
    compproj    = winv(:,compi) * icaEEG(compi,:) - data;    
    squarecomp  = sum(sum(compproj.^2));          
    pvaf(compi) = 100*(1 - squarecomp/squaredata);
end

% Plot the Independent Components
figure
subplot(9, 1, 1:4)
PlotEEG(-data(1:16,:), srate, ChanLabels(1:16), [], 'EEG', 0);
set(gca,'xticklabel',[])
xlabel('')
subplot(9, 1, 5:9)
scale1 = 0.3 * (max(icaEEG(:)) - min(icaEEG(:)));
PlotEEG(-icaEEG(1:20,:), srate, [], scale1, 'Independent Components', 0);
xlabel('')

%%
IC = 1:nbchans;

% Visually identify the independent components representing the artifacts
% to remove.
ICartifacts = [9 12 14 17 20];

IC2use = setdiff(IC, ICartifacts);
PVremoved = sum(pvaf(ICartifacts));

winv = inv(W);
data_ICA = winv(:,IC2use) * icaEEG(IC2use,:);
% This is the cleaned EEG matrix: channels * time points.

% Plot the original and cleaned EEG
time = 580;   % Set the time points to plot
time2plot = ([time time+10]) * srate;   % Plot 10 sec of EEG

figure
subplot(121)
data2plot = data_ICA(1:20,time2plot(1):time2plot(2));
scale = 0.6 * (max(data2plot(:)) - min(data2plot(:)));
PlotEEG(-data(1:20,time2plot(1):time2plot(2)), srate, ChanLabels(1:20), scale, 'EEG', 0);
subplot(122)
PlotEEG(-data2plot, srate, ChanLabels(1:20), scale, 'Cleaned EEG', 0);
