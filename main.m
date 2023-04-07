clear all; close all;
coes = [];
e = earth();
state0 = [7000,0,0,0,sqrt(e.mu/7000),0].';
config_test = struct('dt',10, ...
                     'tspan',[0,6000], ...
                     'state',state0, ...
                     'coes',coes);
sc1 = spacecraft(config_test);
plot3(sc1.state(1,:),sc1.state(2,:),sc1.state(3,:))
hold on 
plot3(state0(1),state0(2),state0(3),"*")
quiver3(state0(1),state0(2),state0(3),state0(4)*1000,state0(5)*1000,state0(6)*1000)
e.plotplanet([0,0,0])
axis equal
xlabel('x'),ylabel('y'),zlabel('z')
grid on
