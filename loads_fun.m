function [A,V,M] = loads_fun(m,cg,I,axial_aero,shear_aero,nx,nz,h,gamma,theta2)
%LOADS_FUN  Calculates internal loads for a rocket structural segment
%   Computes axial force, shear force, and bending moment at the base of
%   a segment due to inertial loads, gravity, and aerodynamic forces.
%
% Inputs:
%   m           : [n, 1] mass of each component in segment [kg]
%   cg          : [n, 1] center of gravity of each component from base [m]
%   I           : [n, 1] moment of inertia of each component [kg⋅m^2]
%   axial_aero  : [1, 1] aerodynamic axial load on segment [N]
%   shear_aero  : [1, 1] aerodynamic shear load on segment [N]
%   nx          : [1, 1] axial load factor (axial acceleration in g's)
%   nz          : [1, 1] lateral load factor (lateral acceleration in g's)
%   h           : [1, 1] height from base to bottom of segment [m]
%   gamma       : [1, 1] flight path angle [rad]
%   theta2      : [1, 1] angular acceleration [rad/s^2]
%
% Outputs:
%   A : [1, 1] axial force at segment base [N]
%   V : [1, 1] shear force at segment base [N]
%   M : [1, 1] bending moment at segment base [N⋅m]
%
% Sign Convention:
%   Shear loads positive toward right
%   Moments positive clockwise

g0 = 9.80665;  % Standard gravity [m/s^2]

% Initialize internal loads
A = 0;  % Axial force
V = 0;  % Shear force
M = 0;  % Bending moment

% Sum contributions from each component in the segment
for i = 1:length(m)
    % Axial force: inertial load + gravity component along flight path
    A = A + m(i)*g0*nx + m(i)*g0*sin(gamma);
    
    % Shear force: lateral inertial + gravity component + rotational effects
    V = V + m(i)*nz*g0 - m(i)*g0*cos(gamma) + m(i)*theta2*(cg(i)-cg(1));
    
    % Apply parallel axis theorem to get moment of inertia about segment base
    I(i) = I(i) + m(i)*(cg(i)-cg(1))^2;
    
    % Bending moment: shear force × moment arm + rotational inertia effects
    M = M + (m(i)*nz*g0 - m(i)*g0*cos(gamma))*(cg(i)-h) ...
          + m(i)*theta2*(cg(i)-cg(1))*(cg(i)-h) ...
          + I(i)*theta2;
end

% Add aerodynamic contributions to total loads
A = A + axial_aero;  % Add aerodynamic drag/thrust
V = V + shear_aero;  % Add aerodynamic lateral force
M = M + shear_aero*(cg(1)-h);  % Add moment from aerodynamic shear

end