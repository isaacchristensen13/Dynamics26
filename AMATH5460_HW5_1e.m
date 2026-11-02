%1d


mu = 0.2;     
rho = 50;     

odetokamuk = @(t, state) [ ...
    state(1)^(0.5) * (state(1)*(1 - mu*state(1)) - state(2)); ...
    rho * (state(1) - 1) * state(2) ...
];


tspan = [0, 80];
x0 = 1.05;
y0 = (1 - mu) * 0.95;


options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);


[t, state] = ode45(odetokamuk, tspan, [x0; y0], options);
x = state(:, 1);
y = state(:, 2);


figure('Position', [100, 100, 900, 650]);


subplot(2, 1, 1);
plot(t, x, 'b-', 'LineWidth', 1.5);
xlabel('Scaled Time \tau');
ylabel('Electron Temperature x(\tau)');
title(sprintf('Sawtooth Oscillation in x(\\tau) (\\mu = %.2f, \\rho = %d)', mu, rho));
grid on;


subplot(2, 1, 2);
plot(x, y, 'r-', 'LineWidth', 1.5); hold on;
plot(1, 1 - mu, 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 7); % Unstable Focus
xlabel('x (Electron Temperature)');
ylabel('y (Instability Amplitude)');
title('Phase Portrait');
legend('Limit Cycle Trajectory', 'Unstable Focus (1, 1-\mu)', 'Location', 'northeast');
grid on; hold off;