
x = [-1:0.1:1]';
n = length(x);
X = [ones(n,1) x x.^2 x.^3];

p = size(X,2);
lambda = 10;

figure
hold on
for i = 1:10
    beta = [randn(p,1)*sqrt(1/lambda)];
    Xbeta = X*beta
    plot(x,Xbeta, 'linewidth',2)
end
xlabel('$x$','fontsize',14, 'interpreter','latex')
ylabel('$E(y) = x^{T}w$','fontsize',14, 'interpreter','latex')
title('10 realizations from 3rd order polynomial with $N(0,\alpha^{-1}I$) prior','fontsize',14, 'interpreter','latex')

print(strcat(figureFolder,'PolynomialPrior'),'-dpng')
