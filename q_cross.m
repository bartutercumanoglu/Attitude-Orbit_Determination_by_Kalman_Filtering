function q_cross = q_cross(q)
    q_cross = [
        q(4)*eye(3)-cross_matrix(q(1:3));
        -q(1:3)'  
    ];
end