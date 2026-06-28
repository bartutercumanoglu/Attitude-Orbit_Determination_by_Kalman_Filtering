function Xdot = f(t,X)

    mu = 398600.4418;   % km^3/s^2
    T_julian = juliandate(datetime(t,'ConvertFrom','posixtime','TimeZone','UTC'));

    r = X(1:3);
    v = X(4:6);

    w = [0; 0; 7.2921159e-5];   % rad/s
    v_rel = v - cross(w,r);

    r_norm = norm(r);

    dr = v;

    RE = 6378.1366;       % km
    J2 = 1.08263e-3;

    x_y_coeff = (3*J2*mu*RE^2/(2*r_norm^5))*(5*(r(3)^2/r_norm^2)-1);
    z_coeff   = (3*J2*mu*RE^2/(2*r_norm^5))*(5*(r(3)^2/r_norm^2)-3);

    a_J2 = [x_y_coeff*r(1);
            x_y_coeff*r(2);
            z_coeff*r(3)];

    a_centralbody = -mu/r_norm^3 * r;

    if norm(v_rel) > 0
        a_drag = -0.5*3.8e-3*(2.2e-7*2.054804/25)*(norm(v_rel)^2)*(v_rel/norm(v_rel));
    else
        a_drag = [0;0;0];
    end

    % Uncomment if you later want 3rd-body
    % r_sun = (planetEphemeris(T_julian,'Earth','Sun'))';
    % r_moon = (planetEphemeris(T_julian,'Earth','Moon'))' - r;
    % a_3rd = (1.9891e30*6.6743e-20/(norm(r_sun)^3))*r_sun + ...
    %         (7.34767309e22*6.6743e-20/(norm(r_moon)^3))*r_moon;

    dv = a_centralbody + a_J2 + a_drag;

    Xdot = [dr; dv];
end