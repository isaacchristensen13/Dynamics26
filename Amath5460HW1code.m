
fig = figure(1);
clf;
set(fig, 'Color', 'w');
ax = axes('Parent', fig);
set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.8 0.8 0.8]);
hold on;


xlim([-2 2]);
ylim([-2 2]);

%Nullclines
f_xnull = @(x,y) -y + x.*(1 - x.^2 - y.^2);
f_ynull = @(x,y)  x + y.*(1 - x.^2 - y.^2);

fimplicit(f_xnull, [-2 2 -2 2], 'r-', 'LineWidth', 2); 
fimplicit(f_ynull, [-2 2 -2 2], 'b-', 'LineWidth', 2); 

%Trajectories
sys = @(t, X) [ -X(2) + X(1)*(1 - X(1)^2 - X(2)^2); ...
                 X(1) + X(2)*(1 - X(1)^2 - X(2)^2) ];

tspan = [0 10];

initial_conditions = [ ...
    0.2,  0.2;   
    0.5, -0.1;   
    1.5,  1.5;   
   -1.8,  0.5;   
   -0.8, -1.5    
];


for i = 1:size(initial_conditions, 1)
    x0 = initial_conditions(i, :);
    [t, X] = ode45(sys, tspan, x0);
    plot(X(:,1), X(:,2), 'k--', 'LineWidth', 1.5);        
    plot(X(1,1), X(1,2), 'ko', 'MarkerFaceColor', 'k');  
end


xlabel('x', 'Color', 'k');
ylabel('y', 'Color', 'k');
title('Phase Space with Nullclines and Trajectories', 'Color', 'k');

leg = legend('x-nullcline (dx/dt=0)', 'y-nullcline (dy/dt=0)', 'Trajectories', 'Location', 'northeast');
set(leg, 'TextColor', 'k', 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);

grid on;
axis square;
hold off;