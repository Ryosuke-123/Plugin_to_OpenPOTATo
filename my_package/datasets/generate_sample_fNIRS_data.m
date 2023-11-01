%% Genarate of sample fNIRS data
% Parameter
Sampling_Freq = 10; % Sampling Frequency[Hz]
Data_length   = 100; % Data Length[s]
CH            = 1;  % Channel Number
t             = 0:1/Sampling_Freq:Data_length - 1/Sampling_Freq;

% Basic fNIRS Signal
baseline = sin(2*pi*0.01*t) + 0.5*sin(2*pi*0.03*t);

% Baseline Shift
baseline_shift_times = sort(rand(1, 5) * T); 
baseline_shift_magnitudes = rand(1, 5) * 0.2; 

for i = 1:length(baseline_shift_times)
    shift_index = round(baseline_shift_times(i) * Sampling_Freq);
    baseline(shift_index:end) = baseline(shift_index:end) + baseline_shift_magnitudes(i);
end

% Spikes
spike_times = sort(rand(1, 10) * Data_length); 
spike_durations = rand(1, 10) * 0.5; 
spike_magnitudes = rand(1, 10) * 0.5;  
for i = 1:length(spike_times)
    spike_index = round(spike_times(i) * Sampling_Freq);
    spike_duration = round(spike_durations(i) * Sampling_Freq);
    spike = spike_magnitudes(i) * sin(2*pi*5*t(spike_index:spike_index+spike_duration)); 
    baseline(spike_index:spike_index+spike_duration) = baseline(spike_index:spike_index+spike_duration) + spike;
end

% Noise fNIRS Signal
Noise_Data = baseline + 1;

%% Glaph
figure(1);
plot(Noise_Data);
