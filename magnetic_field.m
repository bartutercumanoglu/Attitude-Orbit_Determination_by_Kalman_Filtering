function MGM = magnetic_field(X,t)

t_utc = datetime(t,'ConvertFrom','posixtime','TimeZone','UTC');
t_utc_datevec = datevec(t_utc);
r_eci = X(1:3).*1000;
decYear = year(t_utc) + (day(t_utc,'dayofyear') - 1)/365.25;
r_ecef = eci2ecef(t_utc_datevec,r_eci);
r_geo = ecef2lla(r_ecef');
r_geo = r_geo';
B_ned = wrldmagm(r_geo(3),r_geo(1),r_geo(2),decYear,'2025');
phi = deg2rad(r_geo(1));
lam = deg2rad(r_geo(2));
T_ned2ecef = ...
        [-sin(phi)*cos(lam), -sin(lam), -cos(phi)*cos(lam);
         -sin(phi)*sin(lam),  cos(lam), -cos(phi)*sin(lam);
          cos(phi),           0,        -sin(phi)];
B_ecef = T_ned2ecef * B_ned;
B_eci = ecef2eci(t_utc_datevec,B_ecef);
MGM = B_eci;