function s_meas = sun_sensor_model(s_true)
% SUN_SENSOR_MODEL
% ----------------------------------------------------
% Sun sensor measurement model with random angular error
%
% INPUT:
%   s_true : 3x1 true Sun direction vector (unit vector)
%
% OUTPUT:
%   s_meas : 3x1 measured Sun direction vector (unit vector)
%
% Sensor model:
%   - Random angular error
%   - Error angle ∈ [0, 1] degrees
%   - Random rotation axis
%   - Output is normalized
%
% This models a Sun sensor with maximum 1-degree accuracy
% ----------------------------------------------------

    % Ensure unit input
    s_true = s_true / norm(s_true);

    % Maximum error (radians)
    max_err = deg2rad(0.3);

    % Random error angle in [0, max_err]
    delta = max_err * rand();

    % Random rotation axis perpendicular to the Sun direction. This makes
    % the direction error equal to delta instead of smaller by chance.
    axis = randn(3,1);
    axis = axis - dot(axis, s_true) * s_true;
    axis = axis / norm(axis);

    % Rodrigues rotation formula
    K = [   0      -axis(3)  axis(2);
          axis(3)     0     -axis(1);
         -axis(2)  axis(1)     0     ];

    R = eye(3) + sin(delta)*K + (1-cos(delta))*(K^2);

    % Apply rotation
    s_meas = R * s_true;

    % Normalize output
    s_meas = s_meas / norm(s_meas);

end
