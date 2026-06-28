function [q_plus, P_plus, info] = kalman(q_minus, P_minus, q_meas, R, omega, dt, Q)


q_minus = q_minus(:) / norm(q_minus);
q_meas = q_meas(:) / norm(q_meas);

H = eye(3);

q_err = quat_multiply(q_meas, quat_conj(q_minus));
if q_err(4) < 0
    q_err = -q_err;
end

y = 2 * q_err(1:3) / q_err(4);
S = H * P_minus * H' + R;
K = (P_minus * H') / S;

delta_theta = K * y;
P_update = (eye(3) - K * H) * P_minus * (eye(3) - K * H)' + K * R * K';

delta_q = rodrigues_to_quat(delta_theta);
q_update = quat_multiply(delta_q, q_minus);
q_update = q_update / norm(q_update);

P_update = 0.5 * (P_update + P_update');

if nargin >= 7 && ~isempty(omega) && ~isempty(dt) && ~isempty(Q)
    omega = omega(:);
    F = -cross_matrix(omega);
    Phi = eye(3) + F * dt;

    P_plus = Phi * P_update * Phi' + Q * dt;
    P_plus = 0.5 * (P_plus + P_plus');

    q_omega = [omega; 0];
    q_dot = 0.5 * quat_multiply(q_omega, q_update);
    q_plus = q_update + q_dot * dt;
    q_plus = q_plus / norm(q_plus);
else
    q_plus = q_update;
    P_plus = P_update;
end

info.innovation = y;
info.kalman_gain = K;
info.delta_theta = delta_theta;
info.q_error = q_err;

end

function q_conj = quat_conj(q)
q_conj = [-q(1:3); q(4)];
end

function q = rodrigues_to_quat(delta_theta)
q = [delta_theta / 2; 1];
q = q / norm(q);
end

function q = quat_multiply(q_left, q_right)
v_left = q_left(1:3);
s_left = q_left(4);
v_right = q_right(1:3);
s_right = q_right(4);

q = [
    s_left * v_right + s_right * v_left - cross(v_left, v_right);
    s_left * s_right - dot(v_left, v_right)
];
end
