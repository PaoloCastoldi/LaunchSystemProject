function [mv,mnew,cg,CM,Axial_aero,Shear_aero,nz,nx,deltah_ox,deltah_rp1] = aero_plane(const_prop,location,const_aero,mv,r,a,cg,alpha,gamma)
%AERO_PLANE  Calculates aerodynamic loads for unpowered flight (gliding/coasting)
%   Similar to AERODYNAMIC but without thrust vectoring calculations.
%   Used for ballistic coast phases or engine-off conditions.
%
% Inputs:
%   const_prop : [1, 6] propulsion constants [m_dot_lox, m_dot_rp1, m_ox, m_rp1, rho_lox, rho_rp1]
%   location   : [1, 5] flight state [time, altitude, velocity, Mach, dynamic_pressure]
%   const_aero : [1, 6] aerodynamic coefficients [Cn_fin_ant, Cd_fin_ant, Cn_tail, Cd_tail, Cn_nose, Cn_inter]
%   mv         : [n, 1] mass vector of components [kg]
%   r          : [n, 1] radius vector [m]
%   a          : [n, 1] area vector [m^2]
%   cg         : [n, 1] center of gravity positions [m]
%   alpha      : [1, 1] angle of attack [rad]
%   gamma      : [1, 1] flight path angle [rad]
%
% Outputs:
%   mv          : [n, 1] updated mass vector [kg]
%   mnew        : [1, 1] current total vehicle mass [kg]
%   cg          : [n, 1] updated CG positions [m]
%   CM          : [1, 1] overall vehicle center of mass [m]
%   Axial_aero  : [1,11] axial aerodynamic loads per segment [N]
%   Shear_aero  : [1,11] shear aerodynamic loads per segment [N]
%   nz          : [1, 1] lateral load factor
%   nx          : [1, 1] axial load factor
%   deltah_ox   : [1, 1] change in LOX liquid height [m]
%   deltah_rp1  : [1, 1] change in RP-1 liquid height [m]

g0 = 9.80665;  % Standard gravity [m/s^2]

% Extract flight conditions
t = location(1);  % Current time [s]
q = location(5);  % Dynamic pressure [Pa]

% Extract propulsion constants
m_dot_lox_s1 = const_prop(1);   % LOX mass flow rate [kg/s]
m_dot_rp1_s1 = const_prop(2);   % RP-1 mass flow rate [kg/s]
m_ox_stage1 = const_prop(3);    % Initial LOX mass [kg]
m_rp1_stage1 = const_prop(4);   % Initial RP-1 mass [kg]
rho_lox = const_prop(5);        % LOX density [kg/m^3]
rho_rp1 = const_prop(6);        % RP-1 density [kg/m^3]

% Extract aerodynamic coefficients
Cn_fin_ant = const_aero(1);         % Normal force coefficient, forward fins
Cd_fin_ant = const_aero(2);         % Drag coefficient, forward fins
Cn_tail = const_aero(3);            % Normal force coefficient, tail fins
Cd_tail = const_aero(4);            % Drag coefficient, tail fins
Cn_alpha_nose = const_aero(5);      % Normal force slope, nose cone [1/deg]
Cn_alpha_intertage = const_aero(6); % Normal force coefficient, interstage

%% UPDATE PROPELLANT MASS AND CENTER OF GRAVITY
% Calculate remaining propellant (may still be draining even with engines off)
m_ox_stage1 = m_ox_stage1 - m_dot_lox_s1*t;
m_rp1_stage1 = m_rp1_stage1 - m_dot_rp1_s1*t;
m_prop_out = m_dot_lox_s1*t + m_dot_rp1_s1*t;  % Total mass expelled

LOM_effettiva = sum(mv);  % Current mass before update
mv(23) = m_ox_stage1;     % Update LOX mass
mv(26) = m_rp1_stage1;    % Update RP-1 mass

% Calculate change in liquid height
deltah_ox = -m_dot_lox_s1*t/rho_lox/(pi*r(22)^2);
deltah_rp1 = -m_dot_rp1_s1*t/rho_rp1/(pi*r(25)^2);

% Update propellant CG positions
cg_ox_stage1 = cg(23) + deltah_ox;
cg_rp1_stage1 = cg(26) + deltah_rp1;
cg(23) = cg_ox_stage1;
cg(26) = cg_rp1_stage1;

% Calculate overall vehicle center of mass
CM_num = sum(mv.*cg);
mnew = LOM_effettiva - m_prop_out;
CM = CM_num/mnew;

%% CALCULATE AERODYNAMIC FORCES
% Component geometries
D_fwd_skirt = r(4)*2;      % Forward skirt diameter [m]
D_interstage = r(19)*2;    % Interstage diameter [m]

% Normal forces (lift components)
L_nose = Cn_alpha_nose*rad2deg(alpha)*q*pi*(D_fwd_skirt^2)/4;
L_inter = Cn_alpha_intertage*alpha*q*pi*(D_interstage^2-D_fwd_skirt^2)/4;
L_fins_ant = 2*sqrt(2)/2*Cn_fin_ant*q*a(13);  % 4 fins at 45°
L_tails = 2*sqrt(2)/2*Cn_tail*q*a(28);        % 4 fins at 45°

% Drag forces
Drag = [];
for i = [a(1), zeros(1,5), a(19), zeros(1,4)]
    Drag = [Drag, 0.7*q*i];  % 0.7 = empirical correction factor
end
Drag(6) = 2*sqrt(2)/2*Cd_fin_ant*q*a(13);  % Forward fins drag
Drag(11) = 2*sqrt(2)/2*Cd_tail*q*a(28);    % Tail fins drag

%% TRANSFORM TO BODY AXIS
% Convert aerodynamic forces from wind axis to body-fixed axis
Axial_aero = Drag*cos(alpha) + ...
    [-L_nose*sin(alpha), 0, 0, 0, 0, -L_fins_ant*sin(alpha), ...
     -L_inter*sin(alpha), 0, 0, 0, -L_tails*sin(alpha)];

Shear_aero = Drag*sin(alpha) + ...
    [L_nose*cos(alpha), 0, 0, 0, 0, L_fins_ant*cos(alpha), ...
     L_inter*cos(alpha), 0, 0, 0, L_tails*cos(alpha)];

%% CALCULATE LOAD FACTORS (NO THRUST)
% For unpowered flight, load factors from aerodynamics and gravity only

% Lateral load factor: Force balance perpendicular to flight path
shear_loads = sum(Shear_aero) - mnew*g0*cos(gamma);
nz = -shear_loads/(mnew*g0);

% Axial load factor: Force balance along flight path (deceleration only)
nx = (-mnew*g0*sin(gamma) - sum(Axial_aero))/(mnew*g0);

end