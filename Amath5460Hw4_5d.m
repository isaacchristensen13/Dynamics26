



a = 0.2;
b = 0.2;
c = 5.7;


xD0 = 1;
yD0 = 1;
zD0 = 1;

xR0 = 10;
yR0 = -10;
zR0 = 15;


Y0 = [xD0; yD0; zD0; xR0; yR0; zR0];


tspan = [0 50];

options = odeset('RelTol',1e-8,'AbsTol',1e-10);

[t,Y] = ode45(@(t,Y) rossler_system(t,Y,a,b,c), ...
              tspan,Y0,options);


xD = Y(:,1);
yD = Y(:,2);
zD = Y(:,3);

xR = Y(:,4);
yR = Y(:,5);
zR = Y(:,6);




figure;

plot3(xD, yD, zD, 'b-', 'LineWidth', 1.5);
hold on;

plot3(xR, yR, zR, 'r--', 'LineWidth', 1.2);


plot3(xD(1), yD(1), zD(1), 'bo', ...
    'MarkerFaceColor', 'b', 'MarkerSize', 7);

plot3(xR(1), yR(1), zR(1), 'ro', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 7);

xlabel('x');
ylabel('y');
zlabel('z');

title('3D Rössler System Synchronization');
legend('Driver', 'Response', ...
    'Driver initial condition', ...
    'Response initial condition');

grid on;
axis tight;
view(3);
hold off;



function dY = rossler_system(~,Y,a,b,c)

    
    xD = Y(1);
    yD = Y(2);
    zD = Y(3);

    
    xR = Y(4);
    yR = Y(5);
    zR = Y(6);

   
    ex = xR - xD;
    ey = yR - yD;
    ez = zR - zD;

    
    u1 = ey + ez - ex;
    u2 = -ex - (a+1)*ey;
    u3 = -zD*ex - (xD-c)*ez - ex*ez - ez;

    
    dxD = -yD - zD;
    dyD = xD + a*yD;
    dzD = b + zD*(xD-c);

    
    dxR = -yR - zR + u1;
    dyR = xR + a*yR + u2;
    dzR = b + zR*(xR-c) + u3;

    
    dY = [dxD; dyD; dzD; dxR; dyR; dzR];

end

