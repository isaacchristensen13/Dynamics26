
%% 3a



a_values = 0.4:0.01:1.7;


lyapunov_values = zeros(size(a_values));


T = 0.1;            
K = 20000;           
Kburn = 3000;        


initial_condition = [1; 1; 1];


options = odeset('RelTol',1e-9,'AbsTol',1e-11);


for j = 1:length(a_values)

    a = a_values(j);

    
    [t,U] = ode45(@(t,U) system(t,U,a), [0 500], ...
        initial_condition, options);

    idx = find(t >= 400, 1);
    X = U(idx,:)';

    
    v = [1; 0; 0];

    
    total_growth = 0;

    
    for k = 1:(K + Kburn)

        
        initial_augmented = [X; v];

        
        [t,P] = ode45(@(t,P) augmented_system(t,P,a), ...
            [0 T], initial_augmented, options);

        
        final_state = P(end,:)';

        X = final_state(1:3);
        v = final_state(4:6);

        
        growth = norm(v);

        
        v = v/growth;

        
        if k > Kburn
            total_growth = total_growth + log(growth);
        end

    end

    
    lyapunov_values(j) = total_growth/(K*T);

    disp(['a = ', num2str(a), ...
        ', Largest Lyapunov exponent = ', ...
        num2str(lyapunov_values(j))]);

end


figure;
plot(a_values, lyapunov_values, 'b', 'LineWidth', 1.2);
hold on;
yline(0,'k--');
xlabel('a');
ylabel('Largest Lyapunov Exponent \mu_1');
title('Largest Lyapunov Exponent as a Function of a');
grid on;



%% 3b


a = 0.880;


initial_condition = [1; 1; 1];


options = odeset('RelTol',1e-9,'AbsTol',1e-11);


tspan = [0 2700];

[t,U] = ode45(@(t,U) system(t,U,a), tspan, ...
    initial_condition, options);


x = U(:,1);
y = U(:,2);
z = U(:,3);


keep = t >= 700;

x = x(keep);
y = y(keep);
z = z(keep);


figure;
plot3(x,y,z,'b','LineWidth',0.5);
xlabel('x');
ylabel('y');
zlabel('z');
title(['Attractor for a = ', num2str(a)]);
grid on;
axis tight;
view(3);


%% Functions

function dU = system(~,U,a)

x = U(1);
y = U(2);
z = U(3);

dx = y*z;
dy = x-a*y;
dz = 1-x*y;

dU = [dx; dy; dz];

end


function dP = augmented_system(~,P,a)

x = P(1);
y = P(2);
z = P(3);


dx = y*z;
dy = x-a*y;
dz = 1-x*y;


J = [0, z, y;
     1, -a, 0;
    -y, -x, 0];


v = P(4:6);


dv = J*v;


dP = [dx; dy; dz; dv];

end
