%% 2a and 2b

initial_condition = [1; 1; 1];
tspan = [0 200];

[t,U] = ode45(@(t,U) system(t,U), tspan, initial_condition);

x = U(:,1);
y = U(:,2);
z = U(:,3);



keep = t >= 50;

x = x(keep);
y = y(keep);
z = z(keep);

idx = find(t >= 50, 1);

x0 = U(idx,1);
y0 = U(idx,2);
z0 = U(idx,3);


initial_condition1 = [x0; y0; z0];
initial_condition2 = [x0; y0; z0 + 1e-8];


tspan2 = 0:0.01:100;

[t,U1] = ode45(@(t,U) system(t,U), tspan2, initial_condition1);
[t,U2] = ode45(@(t,U) system(t,U), tspan2, initial_condition2);


difference = sqrt(sum((U1-U2).^2,2));


%Attractpr
figure;
plot3(x,y,z,'b','LineWidth',0.5);
xlabel('x');
ylabel('y');
zlabel('z');
title('Numerical Attractor');
grid on;
axis tight;
view(3);

%Log Difference
figure;
plot(t,log(difference),'b','LineWidth',1);
xlabel('Time');
ylabel('ln(||U_1 - U_2||)');
title('Separation of Nearby Trajectories');
grid on;

%% Estimate slope


linear = t >= 20& t <= 90;


p = polyfit(t(linear), log(difference(linear)), 1);


slope = p(1);

disp(['Estimated slope = ', num2str(slope)]);
hold on;
plot(t(linear), polyval(p,t(linear)), 'r', 'LineWidth', 2);
legend('Trajectory separation', 'Linear fit');


%% 2c






T = 0.1;             
K = 5000;            
Kburn = 1000;        


Q = eye(3);


X = [x0; y0; z0];


logR = zeros(K,3);


for k = 1:(K + Kburn)

    
    initial_augmented = [X; Q(:)];

    
    [t,P] = ode45(@(t,P) augmented_system(t,P), ...
        [0 T], initial_augmented);

    
    final_state = P(end,:)';

    X = final_state(1:3);

    
    Pmatrix = reshape(final_state(4:end),3,3);

    
    [Q,R] = qr(Pmatrix);

    
    if k > Kburn
        logR(k-Kburn,:) = log(abs(diag(R)))';
    end

end


lyapunov = sum(logR,1)/(K*T);


lyapunov = sort(lyapunov,'descend');

disp('Estimated Lyapunov spectrum:');
disp(lyapunov);
%% Functions
function dU = system(~,U)

x = U(1);
y = U(2);
z = U(3);

dx = y*z;
dy = x-y;
dz = 1-x*y;

dU = [dx; dy; dz];

end


function dP = augmented_system(~,P)

x = P(1);
y = P(2);
z = P(3);


dx = y*z;
dy = x-y;
dz = 1-x*y;


J = [0, z, y;
    1, -1, 0;
    -y, -x, 0];


Q = reshape(P(4:end),3,3);


dQ = J*Q;


dP = [dx; dy; dz; dQ(:)];

end