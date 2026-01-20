function [alpha,cone_for_table,Ia_n,interstage_for_table,d_D] = CN(location,r,h)
%CN  Calculates angle of attack and aerodynamic parameters for lookup tables
%   Computes angle of attack from crosswind profile and geometric parameters
%   needed for normal force coefficient (CN) table interpolation.
%
% Inputs:
%   location : [1, 5] flight conditions [time, altitude, velocity, Mach, dynamic_pressure]
%   r        : [n, 1] radius of each component [m]
%   h        : [n, 1] height of each component [m]
%
% Outputs:
%   alpha               : [1, 1] angle of attack [rad]
%   cone_for_table      : [1, 1] nose cone parameter for CN lookup
%   Ia_n                : [1, 1] body-to-nose length ratio
%   interstage_for_table: [1, 1] interstage parameter for CN lookup
%   d_D                 : [1, 1] diameter ratio (fwd skirt/interstage)

% Extract flight conditions
altitude = location(2);  % Altitude [km]
v = location(3);         % Velocity [m/s]
Ma = location(4);        % Mach number

% Calculate crosswind velocity based on altitude-dependent wind profile
% Wind profile follows standard atmospheric model for launch vehicle design
if altitude > 35
    vw = 0;  % No significant wind above 35 km
    
elseif altitude < 35 && altitude > 20
    vw = 24.384;  % Constant wind layer [m/s]
    
elseif altitude < 20 && altitude > 14
    % Linear decrease from 76.2 m/s at 14 km to 24.384 m/s at 20 km
    vw = (76.2 - 8.9474*(altitude-14));  % [m/s]
    
elseif altitude < 14 && altitude > 9.6
    vw = 76.2;  % Maximum wind layer [m/s]
    
else
    % Linear increase from ground level
    vw = 6.9288*altitude + 9.144;  % [m/s]
end

% Calculate angle of attack from crosswind
% alpha = arctan(crosswind/flight_velocity)
alpha = atan(vw/v);

% Compute geometric parameters for CN lookup tables
In = h(1);  % Nose cone height [m]
Ia = h(4)+h(7)+h(9)+h(10)+h(12);  % Total body length to interstage [m]
D_fwd_skirt = r(4)*2;  % Forward skirt diameter [m]

% Body-to-nose length ratio (dimensionless)
Ia_n = Ia/In;

% Nose cone parameter: combines compressibility and fineness ratio
% Used for supersonic normal force coefficient lookup
cone_for_table = sqrt(Ma^2-1)*D_fwd_skirt/In;

% Interstage geometry
D_interstage = r(19)*2;  % Interstage diameter [m]
d_D = D_fwd_skirt/D_interstage;  % Diameter ratio

% Interstage cone half-angle
delta_inter = atan((r(19)-r(12))/h(19));

% Interstage parameter: compressibility effect on conical frustum
interstage_for_table = sqrt(Ma^2-1)*sin(delta_inter);

end