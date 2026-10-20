


a = 0.2;
b = 0.2;
c = 5.7;


initial_condition = [1; 1; 1];


tspan = [0 200];


options = odeset('RelTol',1e-9,'AbsTol',1e-11);
[t,U] = ode45(@(t,U) rossler_driver(t,U,a,b,c), ...
    tspan, initial_condition, options);


xD = U(:,1);
yD = U(:,2);
zD = U(:,3);


transient_time = 50;
keep = t >= transient_time;

xD = xD(keep);
yD = yD(keep);
zD = zD(keep);


figure;
plot3(xD,yD,zD,'b','LineWidth',0.5);
xlabel('x_D');
ylabel('y_D');
zlabel('z_D');
title('Rössler System: Driver Chaotic Attractor');
grid on;
axis tight;
view(3);


function dU = rossler_driver(~,U,a,b,c)

xD = U(1);
yD = U(2);
zD = U(3);

dxD = -yD - zD;
dyD = xD + a*yD;
dzD = b + zD*(xD-c);

dU = [dxD; dyD; dzD];

end
