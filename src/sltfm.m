% sltfm
%
% This function calculates the streaming local time-frequency map of a 
% signal
%
% Inputs:
%   signal  : Input signal (time - domain vector)
%   fs      : Sampling frequency of the signal (in Hz)
%   lambda  : Regularization parameter controlling how much a sample 
%             is allowed to deviate from the previous one
%   
% Output:
%   Cs      : Streaming local time-frequency map of the signal

function Cs = sltfm(signal, fs, lambda)
    n_samples = length(signal);    
    if n_samples < 3
        error('sltfm:shortSignal', 'Signal must have at least 3 samples.');        
    end
    %% Adding reflection padding to the signal
    pad_size = n_samples-1;                          
    padded_signal = zeros(1,n_samples + pad_size);
    padded_signal(1,1:pad_size) = flip(signal(1,2:pad_size+1));
    padded_signal(1,pad_size+1:pad_size+n_samples) = signal;
    n_padded = length(padded_signal);
    %%
    % freqs contains the frequencies allowed by the Nyquist criterion
    freqs  = linspace(0,fs/2,n_samples);                                                         
    dt = 1 / fs;            % Sampling period (s)
    n_freqs = length(freqs);
    % Cs represents the streaming local time-frequency map matrix
    Cs = zeros(n_freqs,n_padded);              
    % The time variable is expanded because of the added reflection padding.
    t_padded = 0 : dt : dt * (n_padded-1);    
                                    
    cos_vals = zeros(1,n_padded);   % cos_vals represents cosine matrix
    sin_vals = zeros(1,n_padded);   % sin_vals represents sine matrix
    cos_coef = zeros(1,n_padded);
    sin_coef = zeros(1,n_padded);
    % This loop computes the Fourier coefficients for each frequency 
    % ranging from 0 up to the Nyquist frequency.
    for k = 1:n_freqs   
        % This loop defines the sine and cosine matrices. 
        for idx_t = 1:n_padded          
            cos_vals(1,idx_t) = cos(2*pi*t_padded(1,idx_t)*freqs(1,k));
            sin_vals(1,idx_t) = sin(2*pi*t_padded(1,idx_t)*freqs(1,k));
        end
        % Equation 4-3, Thesis
        cos_coef(1,2) = padded_signal(1,2) / cos_vals(1,2); 
        sin_coef(1,2) = padded_signal(1,2) / sin_vals(1,2); % Equation 4-3, Thesis
        % This loop calculates the streaming local time-frequency map at 
        % frequency 𝑓 = k
        for i = 3:n_padded     
            cos_coef(1,i) = cos_coef(1,i-1) + ...
                (cos_vals(1,i) * (padded_signal(1,i) - cos_vals(1,i)*cos_coef(1,i-1))) / ...
                (cos_vals(1,i)^2 + lambda^2);
            sin_coef(1,i) = sin_coef(1,i-1) + ...
                (sin_vals(1,i) * (padded_signal(1,i) - sin_vals(1,i)*sin_coef(1,i-1))) / ...
                (sin_vals(1,i)^2 + lambda^2);
         
        end
        amp = sqrt((cos_coef .^ 2 + sin_coef .^ 2));
        Cs(k,:) = transpose(amp);
    end
    Cs = Cs(:, pad_size+1 : pad_size+n_samples);
end

