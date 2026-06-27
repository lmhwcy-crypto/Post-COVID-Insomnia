Codes for the Manuscript:

Distinct Electrophysiological Signatures of Post-COVID and Primary Insomnia: Insights from Long-Range Temporal Correlations in EEG

I. ICA.m
1. Load your EEG matrix: 
   record1 = channels * time points
   Set your sampling rate, srate.
   First section is to filter EEG with a Notch filter at 60 Hz.
   record is the filtered EEG: channels * time points.

2. ICA protocol
   record is the notch-filtered EEG: channels * time points.
   Second section is to plot the Independent Components for the EEG matrix.
   Subplot 1 is the original EEG.
   Subplot 2 is the Independent Components.
   
3. Third section is to get the cleaned EEG matrix.
   Visually identify the Independent Components representing the artifacts to remove, and put them into the vector ICartifacts.
   data_ICA is the cleaned EEG matrix: channels * time points.

   Plot the original and cleaned EEG.
   Set the time point to plot and plot 10 seconds of EEG.
   Subplot 1 is the original EEG.
   Subplot 2 is the cleaned EEG.

II. DFA.m
1. Load your EEG matrix: time points * channels * patients.
2. Set your sampling rate, Fs.
3. Set the fit intervals for DFA in seconds, FitInterval. 
   The upper limit should be <10% of signal length.
4. Set the frequency bands of interests, Freqs.
5. Collect the DFA exponent matrix: channels * patients * frequencies.

III. Power_density.m
1. Load your EEG matrix: time points * channels * patients. 
2. Set your sampling rate, Fs.
3. Do Multitaper power spectral density estimate.
4. Set the frequency bands of interest, fidx1, fidx2, fidx3, fidx4, and fidx5.
5. Get the power density matrix, Pwr_COVID, Pwr_Insom, Pwr_NC: Channels * patients * frequencies.