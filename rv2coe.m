function orb_el = rv2coe(X,mu)
I = [1;0;0];
J = [0;1;0];
K = [0;0;1];
R = X(1:3);
V = X(4:6);
r = norm(R);
v = norm(V);
vr = dot(R,V)/r;
H = cross(R,V);
h = norm(H);
N = cross(K,H);
n = norm(N);
i = acos(H(3)/h);
if N(2)>=0
    RAAN = acos(N(1)/n);
else
    RAAN = 2*pi - acos(N(1)/n);
end
E = (1/mu)*(cross(V,H)-mu*(R/r));
e = norm(E);
if E(3)>= 0
    omega = acos(dot(N,E)/(n*e));
else
    omega = 2*pi - acos(dot(N,E)/(n*e));
end

if vr>=0
    theta = acos(dot(E,R)/(e*r));
else
    theta = 2*pi - acos(dot(E,R)/(e*r));
end

orb_el = [h;i;RAAN;e;omega;theta]; 
end

