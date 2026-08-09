% inst_freq
% This function calculates the instantaneous frequency of a signal.
%
% Inputs: 
%   signal  : Input signal (time - domain vector)
%   fs      : Sampling frequency of the signal (in Hz)
%
% Output:
%   w       : Instantaneous frequency of the signal

function w = inst_freq(signal, fs)
    
    n_samples = length(signal);
    if n_samples < 3
        error('inst_freq:shortSignal', ...
            'Signal must have at least 3 samples.');
    end
    
    %% Reflection padding
    pad_size = n_samples - 1;                   
    padded_signal = zeros(1,n_samples + 2*pad_size);
    padded_signal(1,1:pad_size) = flip(signal(1,2:pad_size+1));
    padded_signal(1,pad_size+1:pad_size+n_samples) = signal;
    padded_signal(1,pad_size+n_samples+1:n_samples + 2*pad_size) = ...
        flip(signal(1,n_samples-pad_size:n_samples-1));
    n_padded = length(padded_signal);
    %% 
    analytic_signal = hilbert(padded_signal);             
    hilbert_transform = imag(analytic_signal);                
    % denom = Denominator of equation 3-2, Thesis
    denom = (padded_signal .^2 + hilbert_transform .^2) * 2 * pi;   
    % denom_matrix = Diagonal operator made from the denom
    denom_matrix = diag(denom,0);              
    dt = 1/fs;  % Sampling period (s)
    % d_hilbert = derivative of hilbert_transform
    d_hilbert = diff(hilbert_transform) ./ dt;  
    % d_signal = derivative of padded_signal
    d_signal = diff(padded_signal) ./ dt;
    d_hilbert(1,n_padded) = d_hilbert(1,n_padded-1);
    d_signal(1,n_padded) = d_signal(1,n_padded-1);
    % numer = Numerator in equation 3-2, Thesis
    numer = (padded_signal .* d_hilbert - d_signal .* hilbert_transform);    
    % w = Instantaneous frequency, equation 4-2, Thesis
    w = denom_matrix \ transpose(numer);           
    w = transpose(w(pad_size+1:pad_size+n_samples,1));
end
