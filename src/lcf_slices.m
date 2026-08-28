% lcf_slices
%
% This function calculates the local common-frequency slices of a 2-D
% seismic section at a user-specified frequency.
%
% Inputs:
%   data        : 2-D seismic data (time samples x receivers)
%   fs          : Sampling frequency of the signal (in Hz)
%   win_size    : Smoothing window size (positive integer, should be <<
%                 length(signal))
%   lambda1     : Regularization parameter controlling how much a cosine
%                 coefficient is allowed to deviate from the previous one
%   lambda2     : Regularization parameter controlling how much a sine
%                 coefficient is allowed to deviate from the previous one
%   freq        : Target frequency for which the common-frequency slices is
%                 is computed
%
% Output:
%   cfs         : Local common-frequency slices of the input 2-D data at
%                 the specified frequency, obtained using the given
%                 regularization parameters

function cfs = lcf_slices(data, fs, win_size, lambda1, lambda2, freq)
    dt = 1/fs;                      % Sampling period (s)
    receivers = size(data,2);
    t_sample = length(data(:,1));   % The number of time samples
    t = 0 : dt : dt * (t_sample-1); % Time (s)
    % freqs contains the frequencies allowed by the Nyquist criterion
    freqs  = linspace(0,fs/2,t_sample);   
    n_freqs = length(freqs);   
    C3d = zeros(n_freqs,t_sample,receivers);
    
    %% Creating rectangle smoothing operator
    % smooth_op = Rectangle smoothing operator
    smooth_op = (1/win_size) * eye(t_sample);  
    for i = 1 : win_size-1              
        v = (1/win_size) * ones(1,t_sample-i);
        shift_term = diag (v,-1*i);
        smooth_op = smooth_op + shift_term;
        smooth_op(i,1:i) = 1/i;
    end
    %% Creating triangle smoothing operator
    triangle_op = smooth_op * transpose(smooth_op);                                   
    %% Creating Gaussian smoothing operator, D
    gauss_op = triangle_op * triangle_op;                  
    row_sums = sum(gauss_op, 2);        % Scaling the smoothing operator
    gauss_op = gauss_op ./ row_sums;
    %%
    I = eye(t_sample);  % I = identity operator
    % This loop computes the local common-frequency slice for each 
    % frequency ranging from 0 up to the Nyquist frequency.
    for k = 1:n_freqs
        % cos_diag/sin_diag are the diagonals of the cosine/sine basis 
        % matrices (Psi1/Psi2) at the current frequency
        cos_diag = cos(2*pi*t*freqs(k));
        sin_diag = sin(2*pi*t*freqs(k));
        % Since Psi1/Psi2 are diagonal, transpose(Psi)*Psi reduces to an
        % element-wise square, avoiding a full dense matrix multiplication
        Psi1_sq = diag(cos_diag .^ 2);
        Psi2_sq = diag(sin_diag .^ 2);
        % A1/A2 are the regularized system matrices for the cosine/sine
        % coefficients, combining the smoothing operator (gauss_op) with the
        % lambda-controlled deviation constraint
        A1 = lambda1^2*I + gauss_op*(Psi1_sq - lambda1^2*I);
        A2 = lambda2^2*I + gauss_op*(Psi2_sq - lambda2^2*I);
        % rhs1/rhs2 hold the right-hand side for every receiver at once
        % (data is time_samples x receivers), so a single solve below
        % returns the coefficients for all receivers simultaneously
        rhs1 = gauss_op * (cos_diag' .* data);
        rhs2 = gauss_op * (sin_diag' .* data);
        % a/b = local cosine/sine regression coefficients for all receivers
        % at this frequency (one matrix factorization, reused across
        % receivers since A1/A2 don't depend on the receiver index)
        a = A1 \ rhs1;
        b = A2 \ rhs2;
        % Local amplitude at frequency k for every receiver
        C3d(k,:,:) = sqrt(a.^2 + b.^2);
    end
    cfs = squeeze(C3d(freq,:,:));
end

