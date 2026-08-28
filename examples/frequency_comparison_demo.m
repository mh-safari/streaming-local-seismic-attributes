% frequency_comparision_demo
% This script generates two synthetic signals and compares their
% instantaneous, local, and streaming local frequencies.
%
% Synthetic signals:
%   1. A linear chirp signal
%   2. A signal obtained by convolving a synthetic reflectivity with a 
%      40 Hz Ricker wavelet
%
% The script computes and displays:
%   - Instantaneous frequency (using inst_freq)
%   - Local frequency (using loc_freq)
%   - Streaming local frequency (using stream_lf)
%
% for both synthetic signals.

clc;
clear;
close all;

%% Creating the first synthetic signal (a chirp signal)

fs = 500;                   % Sampling frequency (Hz)
dt = 1 / fs;                % Sampling period (s)
t = 0 : dt : 4;             % Time (s)  
signal1 = chirp(t,5,4,25);  % Creating the chirp signal

%% Creating second synthetic signal (a synthetic seismic trace)

fd = 40;                        % fd = Dominant frequency (Hz)
t_wavelet = -0.2 : dt : 0.2;    %t_wavelet = Time (s)
ricker_wavelet = (1-2*(pi^2)*(fd^2).*(t_wavelet.^2)) ...
    .* exp (-((pi*fd.*t_wavelet).^2));   
reflectivity = zeros(1,140);    % A synthetic reflectivity
reflectivity(1,1) = -1;
reflectivity(1,40) = 0.6;
reflectivity(1,60) = -0.1;
reflectivity(1,80) = -0.5;
reflectivity(1,100) = 0.6;
reflectivity(1,120) = 1;
signal2 = conv(reflectivity,ricker_wavelet);
t2 = 0 : dt : dt*(length(signal2)-1);

%% Instantaneous frequency 
w1 = inst_freq(signal1,fs);
w2 = inst_freq(signal2,fs);

%% Local frequency 
% Window size (wind_size) is chosen larger for signal2 since the seismic 
% trace is more complex than the chirp, requiring stronger smoothing
wl1 = loc_freq(signal1,fs,4);
wl2 = loc_freq(signal2,fs,15);

%% Streaming local frequency
% lambda is chosen larger for signal2 for the same reason: more complex
% or noisier data needs a larger lambda for stronger regularization
ws1 = stream_lf(signal1,fs,0.09);
ws2 = stream_lf(signal2,fs,10.3);

%% Ploting the results 

figure('Name','Two Synthetic Signals') 
  
subplot(4,2,1); plot(t,signal1,LineWidth=1.25)    
title('Chirp Signal', 'FontSize',15)
xlabel('Time(s)','FontSize',15);
ylabel('Amplitude','FontSize',15);
pbaspect([4 1 1])

subplot(4,2,3); plot(t,w1,LineWidth=1.25)     
title('Instantaneous Frequency', 'FontSize',15)
xlabel('Time(s)', 'FontSize',15);
ylabel('Frequency (Hz)', 'FontSize',15);
pbaspect([4 1 1])

subplot(4,2,5); plot(t,wl1,LineWidth=1.25)    
title('Local Frequency', 'FontSize',15)
xlabel('Time(s)', 'FontSize',15);
ylabel('Frequency (Hz)', 'FontSize',15);
pbaspect([4 1 1])

subplot(4,2,7); plot(t,ws1,LineWidth=1.25) 
title('Streaming Local Frequency', 'FontSize',15)
xlabel('Time(s)', 'FontSize',15);
ylabel('Frequency (Hz)', 'FontSize',15);
pbaspect([4 1 1])

subplot(4,2,2); plot(t2,signal2,LineWidth=1.25)    
title('Synthetic Seismic Trace','FontSize',15)
xlabel('Time(s)', 'FontSize',15);
ylabel('Amplitude', 'FontSize',15);
pbaspect([4 1 1])
  
subplot(4,2,4); plot(t2,w2,LineWidth=1.25)     
title('Instantaneous Frequency','FontSize',15)
xlabel('Time(s)','FontSize',15);
ylabel('Frequency (Hz)','FontSize',15);
pbaspect([4 1 1])
  
subplot(4,2,6); plot(t2,wl2,LineWidth=1.25)    
title('Local Frequency','FontSize',15)
xlabel('Time(s)','FontSize',15);
ylabel('Frequency (Hz)','FontSize',15);
axis([0,0.7,15,65]) % zoom in on frequency range of interest
pbaspect([4 1 1])
 
subplot(4,2,8); plot(t2,ws2,LineWidth=1.25)     
title('Streaming Local Frequency','FontSize',15)
xlabel('Time(s)','FontSize',15);
ylabel('Frequency (Hz)','FontSize',15);
axis([0,0.7,15,65]) % zoom in on frequency range of interest
pbaspect([4 1 1])



