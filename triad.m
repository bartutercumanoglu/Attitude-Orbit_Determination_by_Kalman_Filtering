function A = triad(s_i,s_b,b_i,b_b)
t_1_b = s_b;
t_2_b = cross(s_b,b_b);
t_2_b = t_2_b/norm(t_2_b);
t_3_b = cross(t_1_b,t_2_b);

t_1_i = s_i;
t_2_i = cross(s_i,b_i);
t_2_i = t_2_i/norm(t_2_i);
t_3_i = cross(t_1_i,t_2_i);

BT = [t_1_b t_2_b t_3_b];
NT = [t_1_i t_2_i t_3_i];

A = BT*NT';