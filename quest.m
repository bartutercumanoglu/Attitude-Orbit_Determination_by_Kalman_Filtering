function q = quest(ri, bi, a)

B = a(1) * bi(:,1) * ri(:,1)';

for i = 2:size(ri,2)

    B = B + a(i) * bi(:,i) * ri(:,i)';

end

trB = trace(B);

S = B + B';

z = [
    B(2,3)-B(3,2);
    B(3,1)-B(1,3);
    B(1,2)-B(2,1)
];

KB = [
    B+B'-trB*eye(3) z;
    z'             trB
];

[Evec, Eval] = eig(KB);

eigenvalues = diag(Eval);
eigenvectors = Evec;

lamdamax = max(eigenvalues);

rho = lamdamax + trB;

q = [
    adj(rho*eye(3) - S)*z;
    det(rho*eye(3) - S)
];

q = q/norm(q);
