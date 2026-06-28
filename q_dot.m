function q_dot = q_dot(q)
    q_dot = [
        q(4)*eye(3)+cross_matrix(q(1:3));
        -q(1:3)'  
    ];
end