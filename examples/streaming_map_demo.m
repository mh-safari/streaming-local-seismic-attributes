% streaming_map_demo
%
% This script generates a synthetic signal formed by the sum of two
% exponential chirp signals and computes its time-frequency representations
% using three different methods:
%
%   1. Short-Time Fourier Transform (STFT) - using MATLAB's built-in
%   spectrogram function   
%   2. Local Time-Frequency Map - using the ltfm function
%   3. Streaming Local Time-Frequency Map - using sltfm function
%
% The resulting time-frequency map from all three methods are displayed for
% comparison.

clc;
clear;
close all;

dt = 0.004;     % Sampling period (s)
L = 2;          % Record length (s) 2
fs = 1 / dt;    % Sampling frequency (Hz)
t = 0 : dt : L; % Time (s)
n_samples = length(t);  
% freqs contains the frequencies allowed by the Nyquist criterion
freqs  = linspace(0,fs/2,n_samples);               

%% Creating the synthetic signal

signal = chirp(t,10,2,40,"quadratic") + chirp(t,30,2,60,"quadratic");
signal = awgn(signal,1);        % Adding Gaussian noise to the chirp signal

%% STFT map
% Window length (45) and overlap (40) control the time-frequency trade-off:
% a longer window gives better frequency resolution but worse time
% resolution, and vice versa. A large overlap (close to the window length)
% produces a smoother map at the cost of more computation.
stft = abs(spectrogram(signal,45,40));

%% Local time-frequency map 
% win_size (10) sets the smoothing window for the local regression, and
% lambda1/lambda2 (6, 6) set the regularization strength for the cosine/
% sine coefficients. Since Gaussian noise was added to the signal, larger
% values were chosen for both to enforce stronger smoothness and avoid
% overfitting to noise (see Fomel, 2007: noisier data requires larger
% lambda for stronger constraints).
% C represents the local time-frequency map of the signal
C = ltfm(signal,fs,10,6,6);

%% Streaming local time-frequency map
% lambda (5.5) plays the same regularization role as lambda1/lambda2 in
% ltfm, but here it jointly constrains both the cosine and sine streaming
% estimates. It was chosen close to the ltfm lambda values for a fair
% comparison between the batch and streaming versions of the algorithm.
% Cs represents the streaming local time-frequency map of the signal
Cs = sltfm(signal,fs,5.5);  

%% Ploting the results 
figure('Name','Time - Frequency Maps Comparison')
subplot(3,1,1);
imagesc(t,freqs,C);
set(gca,'Ydir','normal')
title('Local Time-Frequency Map','FontSize',15)
xlabel('Time (s)','FontSize',15);
ylabel('Frequency (Hz)','FontSize',15);
colormap("turbo")
clim([0 1])
colorbar
pbaspect([3.5 2.5 1])

subplot(3,1,2);
imagesc(t,freqs,Cs);
set(gca,'Ydir','normal')
title('Streaming Time-Frequency Map','FontSize',15)
xlabel('Time (s)','FontSize',15);
ylabel('Frequency (Hz)','FontSize',15);
colormap("turbo")
clim([0 1])
colorbar
pbaspect([3.5 2.5 1])

subplot(3,1,3);
imagesc(t,freqs,stft); set(gca,'Ydir','normal');
title('STFT','FontSize',15)
xlabel('Time (s)','FontSize',15);
ylabel('Frequency (Hz)','FontSize',15);
colormap("turbo")
colorbar
pbaspect([3.5 2.5 1])