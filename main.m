clear all; close all;
e = earth();
% rp = e.radius+250;
% ecc = 0; i = deg2rad(80); a = rp/(1-ecc);
% RAAN = 0; AOP = 0; TA = 0;


% ISS Example
a = 8059.643886; ecc = 0.17141769; i = deg2rad(27.99996);
RAAN = deg2rad(45.00000175); AOP = deg2rad(30.0127649); TA = deg2rad(39.98711236);

coes = [a,ecc,i,RAAN,AOP,TA].'; % [a,e,i,RAAN,AOP,TA]
% tend = 2*pi/sqrt(e.mu)*(abs(a))^(3/2);
tend = 24*60*60 * 2;
config_test = struct('dt',200, ...
    'tspan',[0,tend], ...
    'state',[], ...
    'coes',coes, ...
    'perts',"j2", ...
    'calc_coes',1, ...
    'drawopt',1);
sc1 = spacecraft(config_test);
state0 = sc1.state(:,1);

figure(1)
plot3(sc1.state(1,:),sc1.state(2,:),sc1.state(3,:))
hold on
plot3(state0(1),state0(2),state0(3),"*")
quiver3(state0(1),state0(2),state0(3),state0(4)*1000,state0(5)*1000,state0(6)*1000)
e.plotplanet([0,0,0])
axis equal
xlabel('x'),ylabel('y'),zlabel('z')
grid on

figure(2)
t = linspace(0,24*60*60 * 2,size(sc1.state,2));
titles = ["a","e","i","\Omega","\omega","\theta"];
ylabels = ["km","e","degrees, \circ","degrees, \circ","degrees, \circ","degrees, \circ"];
for k = 1:size(sc1.coes,1)
    subplot(3,2,k)
    if ismember(k,[3,4,5,6])
        plot(t,rad2deg(sc1.coes(k,:)))
    else
        plot(t,sc1.coes(k,:))
    end
    ylabel(ylabels(k))
    title(titles(k))
    grid on
end
