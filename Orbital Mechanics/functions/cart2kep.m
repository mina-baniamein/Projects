function [semiMajor_axis,ecc,inclination,RAAN,arg_perige,f_0] = cart2kep(r_vec, v_vec, mu, unit)
% this function returns keplerian elements starting from cartesian
% coordinates
% 
% PROTOTYPE:
%   [semiMajor_axis,ecc,inclination,RAAN,arg_perige,f_0] = cart2kep(r_vec, v_vec, mu, unit)
% 
% INPUT: 
%  r_vec[3x1]    cartesian vector of velocity
%  v_vec[3x1]    cartesian vector of velocity
%  
% OUTPUT:
%  semiMajor_axis[1]    semi major axis of the orbit [km]
%  ecc[1]               eccentricity [-]
%  inclination[1]       inclination [rad]
%  RAAN[1]              RAAN (right ascension of the ascending node) [rad]
%  arg_perige[1]        argument fo perigee (right ascension of the ascending node) [rad]
%  f_0[1]               true anomaly
%  unit[str]            'rad' or 'deg'
% 
% 
% CONTRIBUTORS:
%  Pierpaolo Di Carlo 10767871
%  Alessandro Bellezza 10673485
%  Gaia Trovatelli 10582310
%  Mina Baniamein 10627453
%  
% VERSIONS
% Firs version
    r = norm(r_vec); v = norm(v_vec);
    vr = dot(v_vec,r_vec/r); h_vec = cross(r_vec,v_vec);
    h = norm(h_vec);
    semiMajor_axis = -mu*(1/(v^2-2*mu/r));
    e = (1/mu)*cross(v_vec,h_vec)-r_vec/r; ecc=norm(e);
    % i = inclination, inclinazione tra h e z, calcolata come cos^-1(h_z/h)
    i = acos(h_vec(3)/h); %rad
    % se hz>0, i appartiene a [0,pi/2)
    % se hz=0, i = pi/2
    % se hz<0, i appartiene a (pi/2,pi]
    % N = node line, perpendicolare al piano identificato da h e z
    z = [0 0 1]; N = cross(z,h_vec); N_cap = N/norm(N);
    % OMEGA = right ascension of the ascending node, angolo destro tra la
    % gamma-line e l'ascending node
    if N_cap(2) >= 0 && i ~= 0
        OMEGA = acos(N_cap(1)); %rad
    elseif i == 0
        OMEGA = 0;
    else
        OMEGA = 2*pi - acos(N_cap(1)); %rad
    end
    % omega = argument of pericentre, angolo fra N ed e (ovvero il
    % pericentro)
    if e(3) >= 0 && ecc ~= 0
        omega = acos(dot(N_cap,e/ecc)); %rad
    elseif ecc == 0
        omega = 0;
    else
        omega = 2*pi - acos(dot(N_cap,e/ecc)); %rad
    end
    % theta_0 = true anomaly, posizione lungo l'orbita
    if vr >= 0
        theta_0 = acos(dot(r_vec,e)/(r*ecc)); %rad
    else
        theta_0 = 2*pi - acos(dot(r_vec,e)/(r*ecc)); %rad
    end
    switch lower(unit)
        case ('rad')
            RAAN = OMEGA; arg_perige = omega; f_0 = theta_0; inclination = i;
        case ('deg')
            RAAN = rad2deg(OMEGA); arg_perige = rad2deg(omega); f_0 = rad2deg(theta_0); inclination = rad2deg(i);
    end

end