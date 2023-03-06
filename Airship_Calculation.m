clear all ; clc;
 
%% input parameters:

% Operation data
h_service = 500;                    %  service ceiling (meter)
P_L  = 3;                         %  pay load (Kg)
P_h  = 95299;                       %  ambiant pressure at service ceiling (pa)
T_h  = 11.5 + 273;                  %  ambiant temperature at service ceiling  (K)
v    = 20;                          %  expected airship Velocity    (m/s)
AOA  = 6;                           %  AOA of the Airship

% Standard Sea level data
P_SL    = 101325;                      %  ambiant pressure      (pa)
T_Sl    = 288.15;                      %  ambiant pressure      (k)
rho_SL  =1.225;                        %  ambiant air Density   (Kg/m^3)
g       = 9.806;                       %  gravity acceleration  (m/s^2)
mu_Air_SL  = 1.78938e-5;               %  dynamc Viscosity       (N s/m2)
% mu_Air_SL = (0.1458e-5 * T_Sl^(3/2)) / (T_Sl+110.4);     % dynamc Viscosity  (N s/m2)(imperical)

MW_gas = 4;                            % Helium molicular mass  (Kg/Kmol) 

% Geometrical Parameters
a=0.812;            % profile equation parameters
b=8.766;            % profile equation parameters
c=9.997;            % profile equation parameters
d= 6.051;           % profile equation parameters
l=4.0001;          % envelope length                                      

% Systems Specifications
t_material  = 0.002;     % thickness of envelope
prop_eff    = 0.95;     % eff of the prepeller
Pw_control  = 30;       % power of control system  (watt)
t_day       = 12;       % time of day (hour)
t_night     = 12;       % time of night (hour)
eff_bat     = 0.99;      % fuel cell (battery) efficiency 
Pw_sc       = 200;      % assumed power of solar array (watt/m2)
eff_sc      = 0.187;    % fuel cell (battery) efficiency 
As          = 2;        % solar Array Area 


%% calculations theoretical design requirements:

rho_Air     = P_h/(287*T_h);                % Air dinsity at service ceiling    (Kg/m3)
rho_Gas     = P_h/(8314/MW_gas *T_h);       % LTA Gas dinsty at service ceiling (Kg/m3)

V_th        = P_L/(rho_Air - rho_Gas);      % theoretical volume - bouyency theory (m3)
V_Req       = 1.1*V_th;                     % Required volume - with 1.2 safety factor 
M_Gas       = rho_Gas * V_Req;              % required mass of LTA Gas 

P_Gas_SL    = rho_Gas*(8314/MW_gas)*T_Sl;   % Pressure of the Gas inside the ballonates



%% Design Calculatios  (all geometrical calculation considering origin at leading edge)

% Geometrical Design Calculation 
syms x
y = 1/8*sqrt(a*(l-x)*(b*x-l*c^0.5+sqrt(c*l^2-d*l*x)));   % profile curve equation
Dy = diff(y,x);                                          % first derivative of the equation

Ve=double(pi*int(y^2,x,[0,l]));                          % required volume - revolved area theory 
Ae=double(2*pi*int(y*sqrt(1+Dy^2),x,[0,l]));             % required surface area - revolved line theory

x_Cp = real(double(int(x*y,x,[0,l])/int(y,x,[0,l])));    % center of pressure theoreticalcalculation
x_Cg = double(2*pi*int(x*y*sqrt(1+Dy^2),x,[0,l])/Ae);    % center of gravity theoretical calculation
D_Cg = double(subs(y,x_Cg)*2);                           % diameter at center of gravity
 
x_D_max = min(double(solve(Dy==0,x)));                   % max diameter and location
D_max = double(subs(y,x_D_max)*2);
AR = l/D_max;                                            % Aspect ratio of length to max diameter

% hoop stress calculation  at sea level  (am not sure with the results...)
P_static = 0.2308*rho_SL*g*D_Cg*(AR)^2;                         % static pressure on the envelope
P_dynamic = rho_SL* v^2 * Ve*0.44*(sind(2*AOA)/(pi*D_Cg^3));    % dynamic pressure on the envelope
P_diff   =  1.722*g*rho_SL*D_max*AR*sind(AOA);                  % internal differential pressure
P_delta  =  P_static + P_dynamic + P_diff;                      % minimum inner pressure
Hoop_stress = P_delta*D_max/t_material;                         % hoop stresses on teh material to keep it rigid

% Aerodynamic calculations 
Re  =  rho_SL*v*l/mu_Air_SL;                                    % Reynolds number at specified speed at Sea Level
Cdv = (0.172*(AR)^(1/3)+0.252*(AR)^(-1.2)+1.032*(AR)^(-2.7))/(Re)^(1/6);   %Volumetric drag coeff.

Fd=0.5*rho_Air*v^2*Cdv*Ve^(2/3)/2;                              % Drag on the envelope

% power of the propeller
Pw_prop = Fd*v/prop_eff;                                        % power needed to overcome the drag  


% total Power required  
Pw_t = Pw_prop + Pw_control;


% plot the shape of the envelope
fplot(y)
axis([0 l -l l])



