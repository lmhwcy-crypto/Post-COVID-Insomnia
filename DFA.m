% DFA
clear all, close all, clc

% Load your EEG matrix, sleep or waking EEG: time points * channels * patients 
load COVID_notch.mat

nbchans = size(COVID_notch, 2);
nbpts   = size(COVID_notch, 3);
Fs      = 250;     % Your sampling rate
nyquist = Fs/2;

CalcInterval = [0.8 100];         
% Set the fit intervals in seconds. The upper limit should be <10% of
% signal length.
FitInterval  = [2 50];                    
logbin       = 25;
d1  = floor(log10(CalcInterval(1) * Fs));    
d2  = ceil(log10(CalcInterval(2) * Fs));
win = round(logspace(d1, d2, (d2-d1)*logbin));	
winsize = win((CalcInterval(1)*Fs <= win & win <= CalcInterval(2)*Fs));	
DFA_Overlap = 0.5; 

% Frequency bands of intersts: delta, theta, alpha, beta, and gamma
Freqs   = [0.5 4 8 13 30 80];   
nbfreqs = length(Freqs) - 1;

DFA = NaN(length(winsize), nbchans, nbpts, nbfreqs);   
Exp = NaN(nbchans, nbpts, nbfreqs);     % DFA exponents

for fi = 1:nbfreqs
    
    % Filt EEG into the frequency band of interest
    low   = Freqs(fi)
    high  = Freqs(fi + 1);
    order = 1500;
    
    filterweights = fir1(order, [low high]/nyquist);
    
    data_filt  = NaN(size(COVID4_notch));
    data_env   = NaN(size(COVID4_notch));
    
    for pi = 1:nbpts
        data_filt(:,:,pi) = FiltFiltM(filterweights, 1, COVID4_notch(:,:,pi));
        data_env(:,:,pi)  = abs(hilbert(data_filt(:,:,pi)));
    end
    
    for pi = 1:nbpts
        for chi = 1:nbchans
    
        % Start DFA    
        DFA_y = NaN(length(winsize), 1);
        y = data_env(:,chi,pi) ./ mean(data_env(:,chi,pi));  
        y = y - mean(y);
        y = cumsum(y);         		
        for i = 1:length(winsize);				
            D = zeros(floor( length(y)/(winsize(i)*(1-DFA_Overlap)) ), 1);		
            tt = 0;
            for nn = 1 : round(winsize(i)*(1-DFA_Overlap)) : length(y)-winsize(i);	
                tt = tt + 1;
                D(tt) = (mean(fastdetrend(y(nn:nn+winsize(i))).^2, 1)) ^ (1/2);		
            end
            DFA_y(i) = mean(D(1:tt), 1);						
        end  					  	       	
        DFA(:, chi, pi, fi) = DFA_y;
    
        DFA_SmallTimeFit_LogSample = min(find(winsize>=FitInterval(1)*Fs));
        DFA_LargeTimeFit_LogSample = max(find(winsize<=FitInterval(2)*Fs));
        X = [ones(1,DFA_LargeTimeFit_LogSample-DFA_SmallTimeFit_LogSample+1)' log10(winsize(DFA_SmallTimeFit_LogSample:DFA_LargeTimeFit_LogSample))'];
        Y = log10(DFA_y(DFA_SmallTimeFit_LogSample:DFA_LargeTimeFit_LogSample));
        DFA_exp = X\Y;      
        Exp(chi, pi, fi) = DFA_exp(2);    
        % Collect the DFA exponent matrix: channels * patients *
        % frequencies
    
        end
    end
end
