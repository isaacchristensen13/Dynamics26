alpha = 0.5;
N = 200;

x = linspace(-pi, pi, 2000);

W = zeros(size(x));

for n = 0:N
    W = W + 2^(-n*alpha)*cos(2^n*x);
end

plot(x,W)
xlabel('x')
ylabel('W(x)')
title('Approximation of W(x)')
grid on