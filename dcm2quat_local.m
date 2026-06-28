function q = dcm2quat_local(C)

trC = trace(C);

q4 = sqrt((1 + trC)/4);

q3 = (C(1,2) - C(2,1))/(4*q4);

q2 = (C(3,1) - C(1,3))/(4*q4);

q1 = (C(2,3) - C(3,2))/(4*q4);

q = [q1; q2; q3; q4];

end
