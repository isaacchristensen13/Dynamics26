% 4(c)
N = 200;
x = linspace(-1,1,N);

lambda = 0.5;
tol = 1e-8;
maxiter = 100;

initials = {zeros(size(x)), ones(size(x)), cos(pi*x)};

figure
hold on

for j = 1:3

    f = initials{j};
    error = zeros(1,maxiter);

    for n = 1:maxiter

        fnew = zeros(size(x));

        for i = 1:N
            integrand = f ./ (1 + (x(i)-x).^2);
            fnew(i) = sin(2*pi*x(i)) + lambda*trapz(x,integrand);
        end

        error(n) = max(abs(fnew-f));

        if error(n) < tol
            break
        end

        f = fnew;
    end

    semilogy(1:n,error(1:n),'o-')

end

xlabel('Iteration')
ylabel('||f_{n+1}-f_n||_\infty')
title('Fixed Point Iteration for Different Initial States')
legend('f_0 = 0','f_0 = 1','f_0 = cos(\pi x)')
grid on

%%Test different lambda values

lambdas = [0.25 0.5 0.6 0.7 1.0 1.5 2.0 3.0];

for k = 1:length(lambdas)

    lambda = lambdas(k);
    f = zeros(size(x));

    for n = 1:maxiter

        fnew = zeros(size(x));

        for i = 1:N
            integrand = f ./ (1 + (x(i)-x).^2);
            fnew(i) = sin(2*pi*x(i)) + lambda*trapz(x,integrand);
        end

        if max(abs(fnew-f)) < tol
            break
        end

        f = fnew;
    end

    if n < maxiter
        fprintf('lambda = %.1f: converged in %d iterations\n',lambda,n)
    else
        fprintf('lambda = %.1f: did not converge\n',lambda)
    end

end