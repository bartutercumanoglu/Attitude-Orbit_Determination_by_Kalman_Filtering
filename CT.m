  function M = CT(X)
r = X(1:3);
v = X(4:6);
v_r = (dot(r,v)/norm(r))*r/norm(r);
v_p = v-v_r;
i = v_p/norm(v_p);
k = -X(1:3)/norm(X(1:3));
j = cross(k,i);

M_T = [i(1) i(2) i(3);
       j(1) j(2) j(3);
       k(1) k(2) k(3)];
M = M_T';
