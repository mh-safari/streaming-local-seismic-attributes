% loc_freq
% This function calculates the local frequency of a signal.
%
% Inputs: 
%   signal      : Input signal (time - domain vector)
%   fs          : Sampling frequency of the signal (in Hz)
%   win_size    : Smoothing window size (positive integer, should be <<
%                 length(signal))
% Output:
%   wl          : Local frequency of the signal

function wl = loc_freq(signal, fs, win_size)
    
    n_samples = length(signal);
    if n_samples < 3
        error('loc_freq:shortSignal', ...
            'Signal must have at least 3 samples.');
    end
    if win_size < 1 || win_size ~= floor(win_size) || win_size >= n_samples
        error('loc_freq:invalidWindow', ...
            'nb must be a positive integer smaller than the signal length.');
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
    % d_hilbert = derivative of hilbert transform
    d_hilbert = diff(hilbert_transform) ./ dt;  
    % d_signal = derivative of padded_signal   
    d_signal = diff(padded_signal) ./ dt;   
    d_hilbert(1,n_padded) = d_hilbert(1,n_padded-1);
    d_signal(1,n_padded) = d_signal(1,n_padded-1);
    % numer = Numerator in equation 3-2, Thesis
    numer = (padded_signal .* d_hilbert - d_signal .* hilbert_transform);
    
    %% Creating rectangle smoothing operator, smooth_op
    smooth_op = (1/win_size) * eye(n_padded); 
    for i = 1 : win_size-1                           
        v = (1/win_size) * ones(1,n_padded-i);
        shift_term = diag(v,-1*i);
        smooth_op = smooth_op + shift_term;
        smooth_op(i,1:i) = 1/i;
    end
    
    %% Creating triangle smoothing operator, triangle_op
    triangle_op = smooth_op * transpose(smooth_op);                                                   
    row_sums = sum(triangle_op, 2);
    triangle_op = triangle_op ./ row_sums;
    %%
    I = eye(n_padded);  % I = identity operator
    % A natural choice for lambda is the least-squares norm of denom_matrix
    % (Fomel, 2007)
    lambda = norm(denom_matrix,2);
    % wl = Local frequency, equation 19-2, Thesis
    wl = ...
        (lambda^2 * I + triangle_op * (transpose(denom_matrix)*denom_matrix - lambda^2 *I)) ... 
        \ (triangle_op * transpose(denom_matrix) * transpose(numer)); 
    wl = transpose(wl(pad_size+1:pad_size+n_samples,1));

end