clc; clear; close all;

clear magnetometer_model
clear sun_sensor_model
clear kalman

Me = 7.71*10^(12);
We = 7.29*10^(-5);
mu = 398600.4418;
Tilt = deg2rad(9.3);
r0 = [5428.388;2690.028;3271.756];
v0 = [3.70992;0.591109;-6.618145];
X = [r0;v0];

tspan = 1746094405:1:1746094839;

options = odeset('RelTol',1e-12,'AbsTol',1e-12*(ones(1,6)));
[t,X] = ode45(@(t,X) f(t,X),tspan,X,options);                                                                                                                                                                                                                                                      
orb_el = rv2coe(X(:,1),mu);
X = X';
t = t';
Ang = norm(cross(X(1:3,1),X(4:6,1))/norm(X(1:3,1))^2);

T = datetime(t, 'ConvertFrom','posixtime','Format','eeee, d. MMMM uuuu H:mm:ss.SSS');

R = X(1:3,:);
R_norm = vecnorm(R);
alt = R_norm - 6378.1366;

Table = readtable('Ephemeris_data_from_GNSS.xlsx','Sheet','250501');
x = Table.pos_x/1000;
y = Table.pos_y/1000;
z = Table.pos_z/1000;
v_x = Table.vel_x/1000;
v_y = Table.vel_y/1000;
v_z = Table.vel_z/1000;
t_real = Table.Var2;
[t_real_unq, ~, ic]  = unique(t_real,"stable");
[t_real_unq, idx] = sort(t_real_unq);
t_real_unq = t_real_unq';
x_real_u = accumarray(ic,x,[],@mean);
y_real_u = accumarray(ic,y,[],@mean);
z_real_u = accumarray(ic,z,[],@mean);
v_x_real_u = accumarray(ic,v_x,[],@mean);
v_y_real_u = accumarray(ic,v_y,[],@mean);
v_z_real_u = accumarray(ic,v_z,[],@mean);
x_real_u = x_real_u(idx); y_real_u = y_real_u(idx); z_real_u = z_real_u(idx);
v_x_real_u = v_x_real_u(idx); v_y_real_u = v_y_real_u(idx); v_z_real_u = v_z_real_u(idx);
pos_real_u = [x_real_u y_real_u z_real_u]';
v_real_u = [v_x_real_u v_y_real_u v_z_real_u]';
X_real_u = [pos_real_u;v_real_u];
alt_real_u = vecnorm([x_real_u y_real_u z_real_u]');

x_real_interp = interp1(t_real_unq,x_real_u,t,"linear");
y_real_interp = interp1(t_real_unq,y_real_u,t,"linear");
z_real_interp = interp1(t_real_unq,z_real_u,t,"linear");
v_x_real_interp = interp1(t_real_unq,v_x_real_u,t,"linear");
v_y_real_interp = interp1(t_real_unq,v_y_real_u,t,"linear");
v_z_real_interp = interp1(t_real_unq,v_z_real_u,t,"linear");
pos_real_interp = [x_real_interp; y_real_interp; z_real_interp];
v_real_interp = [v_x_real_interp;v_y_real_interp;v_z_real_interp];
alt_real_interp = interp1(t_real_unq,alt_real_u,t,"linear") - 6378.1366;
X_real_interp = [pos_real_interp;v_real_interp];

P_orbit = diag([1 1 1 0.01 0.01 0.01].^2);

R_orbit = diag([ ...
    0.05 0.05 0.05 ...          % km
    0.0005 0.0005 0.0005 ...    % km/s
].^2);

Q_orbit = diag([ ...
    0.001 0.001 0.001 ...
    1e-5 1e-5 1e-5 ...
].^2);

X_orbit_kalman = zeros(6, numel(t));
K_orbit_store = nan(6, 6, numel(t));
innovation_orbit_store = nan(6, numel(t));

P = deg2rad(2)^2 * eye(3);          
R_kalman = deg2rad(0.3)^2 * eye(3);
q_kalman = zeros(4, numel(t));
K_attitude_store = nan(3, 3, numel(t));
innovation_attitude_store = nan(3, numel(t));
Q_kalman = deg2rad(0.06)^2 * eye(3);
P_floor = deg2rad(0.06)^2 * eye(3);
Q_gyro = deg2rad(0.015)^2 * eye(3);
for k = 1:numel(t)
    r = X(1:3,k);
    v = X(4:6,k);
    omega_eci = cross(r, v) / norm(r)^2;   
    orb_el = rv2coe(X(:,k),mu);
    Ang = norm(cross(X(1:3,k),X(4:6,k))/norm(X(1:3,k))^2);
    T_o2i = CT(X(:,k));
    T_julian = juliandate(datetime(t(k),'ConvertFrom','posixtime','TimeZone','UTC'));
    omega_orbit = T_o2i' * omega_eci;
    Sun_Earth = (planetEphemeris(T_julian,'Earth','Sun'))';
    Sun_eci = Sun_Earth - X(1:3,k);
    Sun_orbit = T_o2i' * Sun_eci;
    if k == 1
    X_orbit_kalman(:,k) = X(:,k);
    else
    dt_orbit = t(k) - t(k-1);
    [X_orbit_kalman(:,k), P_orbit, orbit_info(k)] = kalman_orbit( ...
        X(:,k), P_orbit, X_real_interp(:,k), R_orbit, Q_orbit, dt_orbit);
    K_orbit_store(:,:,k) = orbit_info(k).kalman_gain;
    innovation_orbit_store(:,k) = orbit_info(k).innovation;
    end
    Sun_meas = sun_sensor_model(Sun_orbit);
    Sun_eci_store(:,k) = Sun_eci/norm(Sun_eci);
    MGM = (magnetic_field(X(1:3,k),t(k)));
    MGM_orbit = T_o2i' * MGM;
    MGM_store(:,k) = MGM;
    MGM_meas = magnetometer_model(MGM_orbit);
    A = triad(Sun_eci/norm(Sun_eci),Sun_meas/norm(Sun_meas), MGM/norm(MGM),MGM_meas/norm(MGM_meas));
    [angle_1_true, angle_2_true, angle_3_true] = dcm2angle(T_o2i',"ZYX");
    [angle_1_measured, angle_2_measured, angle_3_measured] = dcm2angle(A,"ZYX");
    angles_true_degree(:,k) = rad2deg([angle_1_true; angle_2_true; angle_3_true]);
    angles_measured_degree_triad(:,k) = rad2deg([angle_1_measured; angle_2_measured; angle_3_measured]);
    r_i = [MGM/norm(MGM) Sun_eci/norm(Sun_eci)];
    b_b = [MGM_meas/norm(MGM_meas) Sun_meas/norm(Sun_meas)];
    w = [10 1.0];
    w = w/sum(w);
    q_quest(:,k) = quest(r_i,b_b,w);
    if k == 1
    q_kalman(:,k) = q_quest(:,k);
    else
    dt = t(k) - t(k-1);

    r = X(1:3,k);
    v = X(4:6,k);
    omega_eci = cross(r, v) / norm(r)^2;
    omega_orbit = T_o2i' * omega_eci;

    [q_kalman(:,k), P, kalman_info(k)] = kalman( ...
        q_kalman(:,k-1), P, q_quest(:,k), R_kalman, ...
        omega_orbit, dt, Q_gyro);
    K_attitude_store(:,:,k) = kalman_info(k).kalman_gain;
    innovation_attitude_store(:,k) = kalman_info(k).innovation;
    end
    %q_quest(:,k) = dcm2quat_local(C_bi);
    q_triad(:,k) = dcm2quat_local(A);
    q_real(:,k) = dcm2quat_local(T_o2i');
    C_bi = quat2dcm_local(q_quest(:,k));
    [angle_1_measured, angle_2_measured, angle_3_measured] = dcm2angle(C_bi,"ZYX");
    angles_measured_degree_quest(:,k) = rad2deg([angle_1_measured; angle_2_measured; angle_3_measured]);
end

for k = 1:1746094839-1746094404
    theta(k) = quat_error_angle(q_quest(:,k),q_real(:,k));
end

for k = 1:1746094839-1746094404
    theta1(k) = quat_error_angle(q_triad(:,k),q_real(:,k));
end

for k = 1:numel(t)
    theta_kalman(k) = quat_error_angle(q_kalman(:,k), q_real(:,k));
end

figure
sgtitle("Estimated Position(km) Versus Time(s)")
for k = 1:3
    subplot(3,1,k)
    plot(t-1.746094405000000e+09,R(k,:))
    xlabel('Time (s)');
    if k == 1
        ylabel('X Axis (km)');
    elseif k == 2 
        ylabel('Y Axis (km)');
    else
        ylabel('Z Axis (km)');
    end
    grid on
end

figure
sgtitle('Quaternion Error - Quest (degrees)')
plot(t-1.746094405000000e+09,theta)
xlabel("Time (s)")
ylabel("Angle (degrees)")

figure
sgtitle('Quaternion Error - Triad (degrees)')
plot(t-1.746094405000000e+09,theta1)
xlabel("Time (s)")
ylabel("Angle (degrees)")


figure
sgtitle('Quaternion Error - Kalman Filter - Quest (degrees)')
plot(t-1.746094405000000e+09,theta_kalman)
xlabel("Time (s)")
ylabel("Angle (degrees)")

figure
sgtitle("Attitude Kalman Gain vs Time")
plot(t-1.746094405000000e+09, squeeze(K_attitude_store(1,1,:)), ...
     t-1.746094405000000e+09, squeeze(K_attitude_store(2,2,:)), ...
     t-1.746094405000000e+09, squeeze(K_attitude_store(3,3,:)))
xlabel("Time (s)")
ylabel("Kalman Gain")
legend("K_{\theta x}", "K_{\theta y}", "K_{\theta z}")
grid on

figure
sgtitle("Orbit Kalman Gain vs Time")
subplot(2,1,1)
plot(t-1.746094405000000e+09, squeeze(K_orbit_store(1,1,:)), ...
     t-1.746094405000000e+09, squeeze(K_orbit_store(2,2,:)), ...
     t-1.746094405000000e+09, squeeze(K_orbit_store(3,3,:)))
xlabel("Time (s)")
ylabel("Position Gain")
legend("K_{rx}", "K_{ry}", "K_{rz}")
grid on

subplot(2,1,2)
plot(t-1.746094405000000e+09, squeeze(K_orbit_store(4,4,:)), ...
     t-1.746094405000000e+09, squeeze(K_orbit_store(5,5,:)), ...
     t-1.746094405000000e+09, squeeze(K_orbit_store(6,6,:)))
xlabel("Time (s)")
ylabel("Velocity Gain")
legend("K_{vx}", "K_{vy}", "K_{vz}")
grid on

figure
sgtitle("Attitude Innovation vs Time")
plot(t-1.746094405000000e+09, rad2deg(innovation_attitude_store(1,:)), ...
     t-1.746094405000000e+09, rad2deg(innovation_attitude_store(2,:)), ...
     t-1.746094405000000e+09, rad2deg(innovation_attitude_store(3,:)))
xlabel("Time (s)")
ylabel("Innovation (degrees)")
legend("\theta_x", "\theta_y", "\theta_z")
grid on

figure
sgtitle("Orbit Innovation vs Time")
subplot(2,1,1)
plot(t-1.746094405000000e+09, innovation_orbit_store(1,:), ...
     t-1.746094405000000e+09, innovation_orbit_store(2,:), ...
     t-1.746094405000000e+09, innovation_orbit_store(3,:))
xlabel("Time (s)")
ylabel("Position Innovation (km)")
legend("r_x", "r_y", "r_z")
grid on

subplot(2,1,2)
plot(t-1.746094405000000e+09, innovation_orbit_store(4,:), ...
     t-1.746094405000000e+09, innovation_orbit_store(5,:), ...
     t-1.746094405000000e+09, innovation_orbit_store(6,:))
xlabel("Time (s)")
ylabel("Velocity Innovation (km/s)")
legend("v_x", "v_y", "v_z")
grid on


figure
sgtitle("Position Error (km) vs Time (s)")
for k = 1:3
    subplot(3,1,k)
    plot(t-1.746094405000000e+09,abs(pos_real_interp(k,:)-R(k,:)))
    xlabel("Time (s)")
    if k == 1
        ylabel("X Axis (km)")
    elseif k == 2
        ylabel("Y Axis (km)")
    else
        ylabel("Z Axis (km)")
    end
    grid on
end

figure
sgtitle("Position Error Kalman (km) vs Time (s)")
for k = 1:3
    subplot(3,1,k)
    plot(t-1.746094405000000e+09,abs(X_orbit_kalman(k,:)-R(k,:)))
    xlabel("Time (s)")
    if k == 1
        ylabel("X Axis (km)")
    elseif k == 2
        ylabel("Y Axis (km)")
    else
        ylabel("Z Axis (km)")
    end
    grid on
end


