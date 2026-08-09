% lambda_effect
%
% This script generates a synthetic signal by convolving a Ricker wavelet
% with a synthetic reflectivity. It then computes the streaming local
% frequency of the signal using the stream_lf function for four different
% values of the regularization parameter lambda, and displays the results.

clc;
clear;
close all;

%% Creating a synthetic seismic signal 

fs = 500;                                  % Sampling frequency (Hz)
dt = 1 / fs;                               % Sampling period (s)
fd = 40;                                   % fd = Dominant frequency (Hz)
% t_wavelet = Time axis for the Ricker wavelet, centered at zero 
t_wavelet = -0.2 : dt : 0.2;               
ricker_wavelet = (1-2*(pi^2)*(fd^2).*(t_wavelet.^2)) ...
    .* exp (-((pi*fd.*t_wavelet).^2));   
% reflectivity = A synthetic reflectivity                                                                 
reflectivity = zeros(1,140);               
reflectivity(1,1) = -1;
reflectivity(1,40) = 0.6;
reflectivity(1,60) = -0.1;
reflectivity(1,80) = -0.5;
reflectivity(1,100) = 0.6;
reflectivity(1,120) = 1;
signal = conv(reflectivity,ricker_wavelet); % Creating the synthetic trace
t = 0 : dt : dt*(length(signal)-1);         % t = Time (s)
n_samples = length(signal);   


%% Calculating streaming local frequency
lambda = [0.5, 1, 5, 10];
for j = 1 : 4
    ws = stream_lf(signal, fs, lambda(j));
    subplot(2,2,j);
    plot(t,ws,LineWidth=1.5)
    title(['Streaming local frequency with \lambda = ', num2str(lambda(j))], ...
        'FontSize',14)
    xlabel('Time (s)','FontSize',14);
    ylabel('Frequency (Hz)','FontSize',14);
    axis([0, 0.7, 15, 60])
end

