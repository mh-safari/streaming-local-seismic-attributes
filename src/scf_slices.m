% scf_slices
%
% This function calculates the streaming local common-frequency slices of a
% 2-D seismic section at a user-specified frequency.
%
% Inputs:
%   data        : 2-D seismic data (time samples x receivers)
%   fs          : Sampling frequency of the signal (in Hz)
%   lambda      : Regularization parameter controlling how much a sample
%                 is allowed to deviate from the previous one
%   freq        : Index (1 to n_freqs) of the target frequency for which
%                 the common-frequency slice is computed
%
% Output:
%   scf         : Streaming local common-frequency slice of the input 2-D
%                 data at the specified frequency, obtained using the
%                 given regularization parameter (time samples x receivers)

function scf = scf_slices(data, fs, lambda, freq)

    dt = 1/fs;                       % Sampling period (s)
    t_samples = size(data, 1);       % The number of time samples
    receivers = size(data, 2);       % The number of receivers

    if t_samples < 3
        error('scf_slices:shortSignal', ...
            'Each trace must have at least 3 time samples.');
    end

    % freqs contains the frequencies allowed by the Nyquist criterion
    freqs = linspace(0, fs/2, t_samples);
    n_freqs = length(freqs);

    if freq < 1 || freq > n_freqs || freq ~= floor(freq)
        error('scf_slices:invalidFreq', ...
            'freq must be an integer index between 1 and n_freqs.');
    end

    pad_size = t_samples - 1;
    n_padded = t_samples + pad_size;
    % The time variable is expanded because of the added reflection padding
    t_padded = 0 : dt : dt * (n_padded-1);
    % cos_vals represents cosine matrix
    cos_vals = zeros(n_freqs, n_padded);
    % sin_vals represents sine matrix
    sin_vals = zeros(n_freqs, n_padded); 
    for k = 1:n_freqs
        cos_vals(k,:) = cos(2*pi*t_padded*freqs(k));
        sin_vals(k,:) = sin(2*pi*t_padded*freqs(k));
    end

    Cs3d = zeros(n_freqs, n_padded, receivers);

    for r = 1 : receivers

        %% Adding reflection padding
        signal = transpose(data(:,r));
        padded_signal = zeros(1, n_padded);
        padded_signal(1,1:pad_size) = flip(signal(1,2:pad_size+1));
        padded_signal(1,pad_size+1:pad_size+t_samples) = signal;

        %% Calculating the streaming local time-frequency map
        % Cs represents the streaming local time-frequency map for this
        % receiver
        Cs = zeros(n_freqs, n_padded);
        % cos_coef = streaming cosine coefficients
        cos_coef = zeros(1, n_padded);  
        % sin_coef = streaming sine coefficients
        sin_coef = zeros(1, n_padded);  

        % This loop computes the streaming Fourier coefficients for each
        % frequency ranging from 0 up to the Nyquist frequency.
        for k = 1:n_freqs
            cos_k = cos_vals(k,:);
            sin_k = sin_vals(k,:);

            % Equation 4-3, Thesis
            cos_coef(1,2) = padded_signal(1,2) / cos_k(1,2);
            sin_coef(1,2) = padded_signal(1,2) / sin_k(1,2);

            % This loop calculates the streaming local time-frequency map
            % at frequency 𝑓 = k
            for idx = 3:n_padded
                cos_coef(1,idx) = cos_coef(1,idx-1) + ...
                    (cos_k(1,idx) * (padded_signal(1,idx) - cos_k(1,idx)*cos_coef(1,idx-1))) / ...
                    (cos_k(1,idx)^2 + lambda^2);
                sin_coef(1,idx) = sin_coef(1,idx-1) + ...
                    (sin_k(1,idx) * (padded_signal(1,idx) - sin_k(1,idx)*sin_coef(1,idx-1))) / ...
                    (sin_k(1,idx)^2 + lambda^2);
            end

            amp = sqrt(cos_coef.^2 + sin_coef.^2);
            Cs(k,:) = amp;
        end
        Cs3d(:,:,r) = Cs;
    end
    scf = squeeze(abs(Cs3d(freq, pad_size+1:n_padded, :)));
end
