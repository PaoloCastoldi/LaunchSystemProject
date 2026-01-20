function [T,nz,M_aero] = thrust_eq(shear_aero,cg,cg_gim,CM,m,gamma)
%THRUST_EQ  Calculates thrust required for moment equilibrium (vehicle trim)
%   Determines the thrust magnitude needed to balance aerodynamic moments
%   about the vehicle center of mass during flight.
%
% Inputs:
%   shear_aero : [n, 1] aerodynamic shear loads on each segment [N]
%   cg         : [n, 1] centroid locations from vehicle base [m]
%   cg_gim     : [1, 1] gimbal point location from base [m]
%   CM         : [1, 1] vehicle center of mass location from base [m]
%   m          : [1, 1] total vehicle mass [kg]
%   gamma      : [1, 1] flight path angle [rad]
%
% Outputs:
%   T          : [1, 1] thrust magnitude for equilibrium [N]
%   nz         : [1, 1] normal load factor (lateral g-load)
%   M_aero     : [1, 1] total aerodynamic moment about CM [N⋅m]
%
% Sign Convention:
%   Moments positive clockwise, forces positive toward right

g0 = 9.80665;  % Standard gravity [m/s^2]

% Calculate total aerodynamic moment about center of mass
% Sum moments from distributed shear loads: M = Force × moment_arm
M_aero = 0;
for i = 1:length(cg)
    M_aero = M_aero + shear_aero(i)*(cg(i)-CM);
end

% Solve moment equilibrium: T*(cg_gim-CM) + M_aero = 0
% Thrust at gimbal must counteract aerodynamic moments
eq = @(T) T*(cg_gim-CM) + M_aero;
options = optimoptions('fsolve', 'Display', 'none');
T = fsolve(eq, 1e4, options);

% Calculate normal load factor from force equilibrium
% Total lateral force = aerodynamic + thrust - weight component
shear_loads = sum(shear_aero) + T - m*g0*cos(gamma);
nz = -shear_loads/(m*g0);  % Lateral acceleration in g's

end