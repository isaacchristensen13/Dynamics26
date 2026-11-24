
%% 5a


a = 2.4;
b = 1.4;
c = 0.8;
d = 0.5;


N = 100000;

Nburn = 1000;


x0 = 0.5;
y0 = 0.5;


x = zeros(N+1,1);
y = zeros(N+1,1);

x(1) = x0;
y(1) = y0;


for n = 1:N

    x(n+1) = sin(b*y(n)) + c*sin(b*x(n));
    y(n+1) = sin(a*x(n)) + d*sin(a*y(n));

end


x = x(Nburn+1:end);
y = y(Nburn+1:end);


figure;
plot(x,y,'.','MarkerSize',1);
xlabel('x');
ylabel('y');
title('Attractor of the Two Dimensional Map');
grid on;
axis equal;

%% 5b


epsilon_values = logspace(log10(0.005),log10(0.5),25);


N_values = zeros(size(epsilon_values));


xmin = min(x);
ymin = min(y);


for j = 1:length(epsilon_values)

    epsilon = epsilon_values(j);

    
    box_x = floor((x-xmin)/epsilon);
    box_y = floor((y-ymin)/epsilon);

    
    occupied_boxes = unique([box_x, box_y],'rows');

    N_values(j) = size(occupied_boxes,1);

end


log_inverse_epsilon = log10(1./epsilon_values);
log_N = log10(N_values);


fit_start = 6;
fit_end = 25;


fit_coefficients = polyfit( ...
    log_inverse_epsilon(fit_start:fit_end), ...
    log_N(fit_start:fit_end),1);

box_dimension = fit_coefficients(1);


figure;
plot(log_inverse_epsilon,log_N,'bo','MarkerSize',5);
hold on;


fit_line = polyval(fit_coefficients, ...
    log_inverse_epsilon(fit_start:fit_end));

plot(log_inverse_epsilon(fit_start:fit_end), ...
    fit_line,'r','LineWidth',1.5);

xlabel('log_{10}(1/\epsilon)');
ylabel('log_{10}(N(\epsilon))');
title('Box Counting Dimension');
legend('Calculated values','Linear fit','Location','best');
grid on;


disp(['Estimated box-counting dimension = ', ...
    num2str(box_dimension)]);
