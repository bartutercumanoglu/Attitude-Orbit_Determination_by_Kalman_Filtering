function theta = quat_error_angle(q_est, q_real)
    q_real_conj = [
        -q_real(1:3);
        q_real(4)
    ];
    q_err = norm(q_real)^-2 * [q_cross(q_est) q_est] * q_real_conj;
    q4 = q_err(4);

    theta = 2 * acosd(q4);
end