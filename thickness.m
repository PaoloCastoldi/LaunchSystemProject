function [thick,stress] = thickness(P,M,V,h,r,t,sigma,E,FS,press)
%THICKNESS  Calculates required wall thickness for cylindrical shell
%   Determines minimum thickness to resist axial, bending, shear loads,
%   and buckling. Optionally includes internal pressure effects.
%
% Inputs:
%   P     : [1, 1] axial force [N]
%   M     : [1, 1] bending moment [N⋅m]
%   V     : [1, 1] shear force [N]
%   h     : [1, 1] cylinder height [m]
%   r     : [1, 1] cylinder radius [m]
%   t     : [1, 1] current thickness from design [m]
%   sigma : [1, 1] allowable compressive stress [Pa]
%   E     : [1, 1] Young's modulus [Pa]
%   FS    : [1, 1] safety factor (dimensionless)
%   press : [1, 4] (optional) pressure data [pressure, density, liquid_height, nx]
%
% Outputs:
%   thick  : [1, 3-4] minimum thickness for each failure mode [m]
%   stress : [1, 2-3] stress components [Pa]

tau = 0.5*sigma;  % Shear yield stress (von Mises criterion)
g0 = 9.80665;     % Standard gravity [m/s^2]

if nargin < 10
    %% UNPRESSURIZED CASE: Structure without internal pressure
    
    % Axial stress (compression/tension from axial load)
    sa = P/(2*pi*r*t)*FS;
    
    % Bending stress (from bending moment)
    sb = M*r/(pi*r^3*t)*FS;
    
    % Shear stress (from shear force)
    tshear = V/(2*pi*r*t)*FS;
    
    % Required thickness for combined axial + bending (yield criterion)
    t_ab = (P/(2*pi*r)+(M/(pi*r^2)))/sigma;
    
    % Required thickness for shear
    t_S = V*FS/(2*pi*r*tau);
    
    % Required thickness for buckling (empirical formula for cylindrical shells)
    % Buckling occurs when compressive stress exceeds critical value
    % Function: P/A + M/I - Critical_buckling_stress = 0
    unpress = @(tc) P/(2*pi*r*tc) + M/(pi*r^2*tc) - E*(9*(tc/r)^1.6 + 0.16*(tc/h)^1.3);
    options = optimoptions('fsolve', 'Display', 'none');
    t_c = fsolve(unpress, 0.01, options);
    
    % Output: stresses and required thicknesses
    stress = [abs(sa)+abs(sb), tshear];
    thick = [t_ab, t_S, t_c];
    
else
    %% PRESSURIZED CASE: Structure with internal pressure (propellant tanks)
    
    % Extract pressure data
    p = press(1);       % Internal pressure [Pa]
    rho = press(2);     % Liquid density [kg/m^3]
    hl = press(3);      % Liquid height [m]
    nx = press(4);      % Axial acceleration factor
    
    % Axial stress (compression/tension from axial load)
    sa = P/(2*pi*r*t)*FS;
    
    % Bending stress (from bending moment)
    sb = M*r/(pi*r^3*t)*FS;
    
    % Longitudinal pressure stress (pressure acts to expand cylinder lengthwise)
    spress = (p + rho*g0*nx*hl)*r/(2*t);
    
    % Hoop stress (circumferential stress from internal pressure)
    % This is typically the critical stress for pressure vessels
    shoop = (p + rho*g0*nx*hl)*r/(t);
    
    % Shear stress (from shear force)
    tshear = V/(2*pi*r*t)*FS;
    
    % Required thickness for combined axial + bending - pressure relief
    % Pressure reduces compressive stress requirement
    t_ab = FS/sigma*(P/(2*pi*r) + (M/(pi*r^2)) - r/2*(p + rho*g0*nx*hl))/sigma;
    
    % Required thickness for shear
    t_S = V*FS/(2*pi*r*tau);
    
    % Required thickness for buckling with pressure stiffening
    % Internal pressure increases buckling resistance
    press_eq = @(tc) P/(2*pi*r*tc) + M/(pi*r^2*tc) - ...
        (9*(tc/r)^0.6 + 0.16*(tc/r)^0.3*(r/h)^1.3 + min([0.191*p/E*(r/tc)^2, 0.229]))*E*tc/r;
    options = optimoptions('fsolve', 'Display', 'none');
    t_c = fsolve(press_eq, 0.005, options);
    
    % Required thickness for hoop stress (primary pressure containment)
    t_ho = FS/sigma*r*(p + rho*g0*nx*hl);
    
    % Output: stresses and required thicknesses
    stress = [abs(sa)+abs(sb)-spress, shoop, tshear];
    thick = [t_ab, t_S, t_ho, t_c];
end

end