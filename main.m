clear all; close all;

addpath("planets")
e = earth();
% % Hyperbolic Orbit with Inclination Example
% rp = e.radius+500;
% ecc = 1.2; i = deg2rad(25); a = rp/(1-ecc);
% RAAN = 0; AOP = 0; TA = 0;

% % Circular Orbit Example
% a = e.radius + 100; ecc = 0; i = deg2rad(0); 
% RAAN = 0; AOP = 0; TA = 0;

% ISS Example
a = 8059.643886; ecc = 0.17141769; i = deg2rad(27.99996);
RAAN = deg2rad(45.00000175); AOP = deg2rad(30.0127649); TA = deg2rad(39.98711236);

coes = [a,ecc,i,RAAN,AOP,TA].'; % [a,e,i,RAAN,AOP,TA]
% tend = 2*pi/sqrt(e.mu)*(abs(a))^(3/2) * 10; % period time * number of orbits
tend = 24*60*60 * 3; % 2 days
config_test = struct('dt',100, ...
    'tspan',[0,tend], ...
    'state',[], ...
    'coes',coes, ...
    'perts',["j2","aero"], ...
    'calc_coes',1, ...
    'animateopt',1, ...
    'frameskip',50, ...
    'solver','rk45');
specs_test = struct();
sc1 = spacecraft(config_test,specs_test); % cool stuff done here
