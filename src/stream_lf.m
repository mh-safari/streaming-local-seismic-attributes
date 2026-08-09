% stream_lf
% This function calculates the streaming local frequency of a signal.
%
% Inputs: 
%   signal      : Input signal (time - domain vector)
%   fs          : Sampling frequency of the signal (in Hz)
%   lambda      : Regularization parameter controlling how much a sample 
%                 is allowed to deviate from the previous one
% Output:
%   ws          : Streaming local frequency of the signal

function ws = stream_lf(signal, fs, lambda)
    n_samples = length(signal);
    if n_samples < 3
        error('stream_lf:shortSignal', 'Signal must have at least 3 samples.');
    end
    %% Reflection padding, rp = padding size
    
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
    dt = 1/fs;  % Sampling period (s)
    % d_hilbert = derivative of hilbert_transform
    d_hilbert = diff(hilbert_transform) ./ dt; 
    % d_signal = derivative of padded_signal
    d_signal = diff(padded_signal) ./ dt;
    d_hilbert(1,n_padded) = d_hilbert(1,n_padded-1);
    d_signal(1,n_padded) = d_signal(1,n_padded-1);
    % numer = Numerator in equation 3-2, Thesis
    numer = (padded_signal .* d_hilbert - d_signal .* hilbert_transform);                   
    
    %% Equation 4-3, Thesis
    ws = zeros(1,n_padded);   % ws = Streaming local frequency
    ws(1,1) = numer(1,1) / denom(1,1);                
    
    for i = 2:n_padded
        ws(1,i) = ws(1,i-1) + ...
            (denom(1,i) * (numer(1,i) - denom(1,i)*ws(1,i-1))) / ...
            (denom(1,i)^2 + lambda^2);
    end
    ws = ws(1,pad_size+1:pad_size+n_samples);
    %%
end