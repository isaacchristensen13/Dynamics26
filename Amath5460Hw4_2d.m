



alpha = 1;
beta = 0.5;

% Critical value
epsilon0 = 2*alpha*beta^(3/2);

% Equilibrium
xstar = acos(beta);
ystar = 0;
zstar = 1/sqrt(beta);
equilibrium = [xstar; ystar; zstar];




epsilon_values = [0.8, 0.65];


initial_condition = equilibrium + [0.05; 0.01; 0.02];


tspan = [0 400];

options = odeset('RelTol',1e-8,'AbsTol',1e-10);

figure;

for k = 1:2
    epsilon = epsilon_values(k);

    
    [t,U] = ode45(@(t,U) governor(t,U,alpha,beta,epsilon), ...
                  tspan, initial_condition, options);

    x = U(:,1);
    y = U(:,2);
    z = U(:,3);

    
    subplot(2,2,k);
    plot(t,x,'LineWidth',1.2);
    hold on;
    yline(xstar,'--');
    xlabel('Time');
    ylabel('x(t)');
    title(sprintf('\\epsilon = %.2f',epsilon));
    grid on;

    
    subplot(2,2,k+2);
    plot(x,y,'LineWidth',1.2);
    hold on;
    plot(xstar,ystar,'ro','MarkerFaceColor','r');
    xlabel('x');
    ylabel('y');
    title(sprintf('Phase portrait, \\epsilon = %.2f',epsilon));
    grid on;
end

function dU = governor(~,U,alpha,beta,epsilon)
    x = U(1);
    y = U(2);
    z = U(3);

    dx = y;
    dy = sin(x)*(z^2*cos(x)-1) - epsilon*y;
    dz = alpha*(cos(x)-beta);

    dU = [dx; dy; dz];
end
