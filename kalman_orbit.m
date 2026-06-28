function [x_plus, P_plus, info] = kalman_orbit(x_pred, P_minus, z_gnss, R, Q, dt)


x_pred = x_pred(:);
z_gnss = z_gnss(:);

if numel(x_pred) ~= 6 || numel(z_gnss) ~= 6
    error('kalman_orbit:InputSize', 'x_pred and z_gnss must both be 6x1 states.');
end

if nargin >= 6 && ~isempty(Q) && ~isempty(dt)
    P_minus = P_minus + Q * dt;
elseif nargin >= 5 && ~isempty(Q)
    P_minus = P_minus + Q;
end

H = eye(6);
y = z_gnss - H * x_pred;
S = H * P_minus * H' + R;
K = (P_minus * H') / S;

x_plus = x_pred + K * y;
P_plus = (eye(6) - K * H) * P_minus * (eye(6) - K * H)' + K * R * K';
P_plus = 0.5 * (P_plus + P_plus');

info.innovation = y;
info.kalman_gain = K;
info.residual_norm_position = norm(y(1:3));
info.residual_norm_velocity = norm(y(4:6));

end
