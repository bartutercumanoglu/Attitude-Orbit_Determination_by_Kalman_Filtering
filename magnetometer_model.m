function b_meas = magnetometer_model(b_true)
% MAGNETOMETER_MODEL
% ----------------------------------------------------
% Magnetometer measurement model with random error
%
% INPUT:
%   b_true : 3x1 true magnetic field vector in BODY frame
%            (unit or non-unit)
%
% OUTPUT:
%   b_meas : 3x1 measured magnetic field vector (unit)
%
% Error model:
%   - Random bias (slowly varying / constant per call)
%   - Random Gaussian noise
%   - Output normalized for attitude algorithms
% ----------------------------------------------------

    % Ensure column vector
    b_true = b_true(:);

    % --- USER-DEFINED SENSOR PARAMETERS ---
    bias_max   = 100;    % maximum constant bias magnitude (nT equivalent)
    noise_std  = 20;     % 1-sigma white noise (nT equivalent)

    % --- CONSTANT SENSOR BIAS ---
    % A real magnetometer bias changes slowly, so keep it fixed between calls.
    persistent bias
    if isempty(bias)
        bias_dir = randn(3,1);
        bias_dir = bias_dir / norm(bias_dir);
        bias_mag = bias_max * rand();      % [0, bias_max]
        bias     = bias_mag * bias_dir;
    end

    % --- RANDOM NOISE ---
    noise = noise_std * randn(3,1);

    % --- APPLY ERROR ---
    b_meas = b_true + bias + noise;

    % --- NORMALIZE OUTPUT ---
    b_meas = b_meas / norm(b_meas);

end
