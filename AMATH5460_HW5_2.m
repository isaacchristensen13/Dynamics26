%% 2a
a = 0.2;
b = 0.2;
c = 5.7;

initial_condition = [1; 1; 1];
tspan = [0 200];

[t,U] = ode45(@(t,U) rossler(t,U,a,b,c), tspan, initial_condition);

x = U(:,1);
y = U(:,2);
z = U(:,3);


keep = t >= 50;

x = x(keep);
y = y(keep);
z = z(keep);


figure;
plot3(x,y,z,'b','LineWidth',0.5);
xlabel('x');
ylabel('y');
zlabel('z');
title('Rössler Attractor');
grid on;
axis tight;
view(3);




%% 2b



a = 0.2;
b = 0.2;
c = 5.7;

initial_condition = [1; 1; 1];
tspan = [0 2000];

options = odeset('MaxStep',0.05);

[t,U] = ode45(@(t,U) rossler(t,U,a,b,c), ...
    tspan, initial_condition, options);

x = U(:,1);
y = U(:,2);
z = U(:,3);


transient_time = 50;

x_cross = [];
z_cross = [];

for i = 1:length(t)-1

    
    if y(i) > 0 && y(i+1) <= 0

        
        alpha = y(i)/(y(i)-y(i+1));

        t_intersect = t(i) + alpha*(t(i+1)-t(i));
        x_intersect = x(i) + alpha*(x(i+1)-x(i));
        z_intersect = z(i) + alpha*(z(i+1)-z(i));

        
        if t_intersect >= transient_time && x_intersect < 0

            x_cross = [x_cross; x_intersect];
            z_cross = [z_cross; z_intersect];

        end
    end
end

NTest = length(x_cross);

fprintf('Number of intersections collected: %d\n', NTest);

if NTest < 200
    error('Not enough intersections. Increase simulation time.');
end
N = 200;

x_cross = x_cross(1:N);
z_cross = z_cross(1:N);



%% 2c




figure;
plot(x_cross, z_cross, 'b.', 'MarkerSize', 10);
xlabel('x_n');
ylabel('z_n');
title('Poincare Section: z_n vs. x_n');
grid on;


figure;
plot(x_cross(1:N-1), x_cross(2:N), 'r.', 'MarkerSize', 10);
xlabel('x_n');
ylabel('x_{n+1}');
title('Poincare Map: x_{n+1} vs. x_n');
grid on;

%% 2d



degree = 10;


[Q_coeff,~,muQ] = polyfit(x_cross, z_cross, degree);


[ P_coeff,~,muP] = polyfit(x_cross(1:N-1), ...
    x_cross(2:N), degree);


x_fit = linspace(min(x_cross), max(x_cross), 500);

figure;
plot(x_cross, z_cross, 'b.', 'MarkerSize', 10);
hold on;
plot(x_fit, polyval(Q_coeff,x_fit,[],muQ), 'r', 'LineWidth', 1.5);
xlabel('x_n');
ylabel('z_n');
title('Polynomial Approximation Q(x)');
legend('Data','Polynomial Fit');
grid on;


figure;
plot(x_cross(1:N-1), x_cross(2:N), 'b.', 'MarkerSize', 10);
hold on;
plot(x_fit, polyval(P_coeff,x_fit,[],muP), 'r', 'LineWidth', 1.5);
xlabel('x_n');
ylabel('x_{n+1}');
title('Polynomial Approximation P(x)');
legend('Data','Polynomial Fit');
grid on;

%% 2e




xmin = min(x_cross);
xmax = max(x_cross);


x_grid = linspace(xmin, xmax, 10000);

for k = 1:3

    
    f = @(x) iterate_P(x,k,P_coeff,muP,xmin,xmax) - x;

    roots_found = [];

    
    for i = 1:length(x_grid)-1

        f1 = f(x_grid(i));
        f2 = f(x_grid(i+1));

        if isfinite(f1) && isfinite(f2)

            if f1*f2 < 0
                root = fzero(f,[x_grid(i),x_grid(i+1)]);
                roots_found = [roots_found; root];
            end

        end
    end

    
    roots_found = uniquetol(roots_found,1e-6);

   
    valid_roots = [];

    for j = 1:length(roots_found)

        x0 = roots_found(j);

        if k == 1
            valid_roots = [valid_roots; x0];

        elseif k == 2
            if abs(iterate_P(x0,1,P_coeff,muP,xmin,xmax)-x0) > 1e-4
                valid_roots = [valid_roots; x0];
            end

        elseif k == 3
            if abs(iterate_P(x0,1,P_coeff,muP,xmin,xmax)-x0) > 1e-4 && ...
                    abs(iterate_P(x0,2,P_coeff,muP,xmin,xmax)-x0) > 1e-4
                valid_roots = [valid_roots; x0];
            end
        end
    end

    
    fprintf('\nPeriod %d:\n',k);

    if isempty(valid_roots)
        fprintf('No periodic points found.\n');
    end

    for j = 1:length(valid_roots)

        x0 = valid_roots(j);
        z0 = polyval(Q_coeff,x0,[],muQ);

        fprintf('x = %.8f, z = %.8f\n',x0,z0);
        fprintf('Initial condition: (%.8f, 0, %.8f)\n',x0,z0);

    end

end

%% 2e confirm periodic orbits

test_roots = [-7.21143453;
    -8.19356138;
    -4.65849544;
    -8.94454140;
    -8.63850359;
    -6.47806210;
    -5.35305963;
    -3.56693852;
    -2.92431266];

fprintf('\nPoincare map verification:\n');

for i = 1:length(test_roots)

    x0 = test_roots(i);
    x_next = polyval(P_coeff,x0,[],muP);

    fprintf('P(%.6f) = %.6f\n',x0,x_next);

end




%% 2f




x_initial = [-7.21143453, -8.19356138, ...
    -8.94454140, -8.63850359];

orbit_names = {'Period 1 Orbit', 'Period 2 Orbit', ...
    'Period 3 Orbit 1', 'Period 3 Orbit 2'};


tspan = [0 20];
options = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',0.01);

for i = 1:length(x_initial)

    
    z_initial = polyval(Q_coeff,x_initial(i),[],muQ);

    
    initial_condition = [x_initial(i); 0; z_initial];

    
    [t,U] = ode45(@(t,U) rossler(t,U,a,b,c), ...
        tspan, initial_condition, options);

    
    figure;
    plot3(U(:,1),U(:,2),U(:,3),'b','LineWidth',1);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title(orbit_names{i});
    grid on;
    axis tight;
    view(3);

end
%% Functions
function dU = rossler(~,U,a,b,c)

x = U(1);
y = U(2);
z = U(3);

dx = -y - z;
dy = x + a*y;
dz = b + z*(x-c);

dU = [dx; dy; dz];

end

function result = iterate_P(x,k,P_coeff,muP,xmin,xmax)

result = x;

for i = 1:k

    result = polyval(P_coeff,result,[],muP);


    if result < xmin || result > xmax
        result = NaN;
        return;
    end

end

end