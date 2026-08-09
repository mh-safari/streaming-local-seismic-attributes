% ltfm
%
% This function calculates the local time-frequency map of a signal
%
% Inputs:
%   signal      : Input signal (time - domain vector)
%   fs          : Sampling frequency of the signal (in Hz)
%   win_size    : Smoothing window size (positive integer, should be <<
%                 length(signal))
%   lambda1     : Regularization parameter controlling how much a cosine 
%                 sample is allowed to deviate from the previous one              
%   lambda2     : Regularization parameter controlling how much a sine 
%                 sample is allowed to deviate from the previous one
% Output:
%   C           : local time-frequency map of the signal

function C = ltfm(signal, fs, win_size, lambda1, lambda2)

    n_samples = length(signal);
    if n_samples < 3
        error('ltfm:shortSignal', 'Signal must have at least 3 samples.');
    end
    if win_size < 1 || win_size ~= floor(win_size) || win_size >= n_samples
        error('ltfm:invalidWindow', ...
            'nb must be a positive integer smaller than the signal length.');
    end
    % fre contains the frequencies allowed by the Nyquist criterion
    freqs  = linspace(0,fs/2,n_samples);           
    dt = 1 / fs;                % Sampling period (s)
    t = 0 : dt : dt*(n_samples-1); % Time (s)
    n_freqs = length(freqs);                                          
    psi1 = zeros(n_freqs,n_freqs);          % psi1 represents cosine matrix
    psi2 = zeros(n_freqs,n_freqs);          % psi2 represents sine matrix
    % C represents the local time-frequency map matrix
    C = zeros(n_freqs,n_samples);   
 
    %% Creating rectangle smoothing operator
    % smooth_op = Rectangle smoothing operator
    smooth_op = (1/win_size) * eye(n_samples);    
    for i = 1 : win_size-1            
        v = (1/win_size) * ones(1,n_samples-i);
        shif_term = diag(v,-1*i);
        smooth_op = smooth_op + shif_term;
        smooth_op(i,1:i) = 1/i;
    end
   
    %% Creating triangle smoothing operator
    triangle_op = smooth_op * transpose(smooth_op);                                   
    
    %% Creating Gaussian smoothing operator, D
    gauss_op = triangle_op * triangle_op;                  
    row_sums = sum(gauss_op, 2);        % Scaling the smoothing operator
    gauss_op = gauss_op ./ row_sums;
    %%
    I = eye(n_samples);            % I = identity operator
    % This loop computes the Fourier coefficients for each frequency 
    % ranging from 0 up to the Nyquist frequency.
    for k = 1:n_freqs                 
        % This loop defines the sine and cosine  matrices. 
        for time = 1:n_samples                   
            psi1(time,time) = cos(2*pi*t(1,time)*freqs(1,k));
            psi2(time,time) = sin(2*pi*t(1,time)*freqs(1,k));
        end
        % These steps calculate the local time-frequency map at frequency 𝑓 = 𝑘
        a = (lambda1^2 * I + gauss_op * (transpose(psi1)*psi1- lambda1^2 *I)) \ ...
            (gauss_op * transpose(psi1) * transpose(signal));
        b = (lambda2^2 * I + gauss_op * (transpose(psi2)*psi2- lambda2^2 *I)) \ ...
            (gauss_op * transpose(psi2) * transpose(signal));
        c = sqrt((a .^ 2 + b .^ 2));
        C(k,:) = transpose(c);
    end
end