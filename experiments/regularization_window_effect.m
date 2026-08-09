% regularization_window_effect
%
% This script generates a synthetic signal by convolving a Ricker wavelet
% with a synthetic reflectivity. It then computes the local frequency of
% the signal using the loc_freq function for four different values of the
% regularization window size, and displays the results.

clc;
clear;
close all;

%% Creating a synthetic seismic signal

fs = 500;                                   % Sampling frequency (Hz)
dt = 1 / fs;                                % Sampling period (s)
fd = 40;                                    % fd = Dominant frequency (Hz)
% t_wavelet = Time axis for the Ricker wavelet, centered at zero 
t_wavelet = -0.2 : dt : 0.2;        
ricker_wavelet = (1-2*(pi^2)*(fd^2).*(t_wavelet.^2)) ...
    .* exp (-((pi*fd.*t_wavelet).^2)); 
reflectivity = zeros(1,140);                      
reflectivity(1,1) = -1;
reflectivity(1,40) = 0.6;
reflectivity(1,60) = -0.1;
reflectivity(1,80) = -0.5;
reflectivity(1,100) = 0.6;
reflectivity(1,120) = 1;
signal = conv(reflectivity,ricker_wavelet);  % Creating the synthetic trace
n_samples = length(signal);         
t = 0 : dt : dt*(length(signal)-1);         % t = Time (s)

%% Calculating local frequency
% Spans from minimal to strong smoothing
win_size = [2, 5, 20, 50];  
for j = 1:4
    wl = loc_freq(signal, fs, win_size(j)); % wl = Local frequency
    subplot(2,2,j);
    plot(t,wl,LineWidth=1.5)
    title(['Local frequency with window size = ', num2str(win_size(j)) ,' samples'], ...
        'FontSize',14)
    xlabel('Time (s)','FontSize',14);
    ylabel('Frequency (Hz)','FontSize',14);
    axis([0, 0.7, 30, 50])
end