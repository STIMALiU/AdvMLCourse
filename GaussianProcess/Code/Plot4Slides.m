%% Generate plots for the lecture slides for the GP part of AdvML course
% Author: Mattias Villani, http://mattiasvillani.com

%% Input
set(0,'DefaultTextInterpreter', 'latex')
nGridPoints = 100;
jitter = 10^-100;  % Adding little numbers to the diagonal of the covariance matrix for stability
figureFolder = '/Users/mv/Dropbox/Teaching/AdvML/GaussianProcess/Slides/Graphics/';
fontSize = 12;
symbols = {'b','r','g'}; % Plot Symbols

%% Bivariate normal
subplot(2,2,1) 
sigmas = [1 1];
rho = 0.0;
PlotBivariateNormalWithMargins([0 0], sigmas, rho, 'off', 0)
title('$\rho=0, \sigma_1 = 1, \sigma_2 = 1$')

subplot(2,2,2) 
sigmas = [1 1];
rho = 0.5;
PlotBivariateNormalWithMargins([0 0], sigmas, rho, 'off', 0)
title('$\rho=0.5, \sigma_1 = 1, \sigma_2 = 1$')


subplot(2,2,3) 
sigmas = [1 1];
rho = 0.95;
PlotBivariateNormalWithMargins([0 0], sigmas, rho, 'off', 0)
title('$\rho=0.95, \sigma_1 = 1, \sigma_2 = 1$')

subplot(2,2,4) 
sigmas = [0.25 1];
rho = 0.5;
PlotBivariateNormalWithMargins([0 0], sigmas, rho, 'off', 0)
title('$\rho=0.5, \sigma_1 = 1/4, \sigma_2 = 1$')


print(strcat(figureFolder,'BivariateNormalContours'),'-dpng')


%% Nonparametric = one parameter at every point

xGrid = linspace(-1,1,nGridPoints)';
figure('name','One parameter at every point')
hold on

% Square exp
K = KernelFunctions('SquaredExp', xGrid, xGrid, [1 0.1]);

fSim = mvnrnd(zeros(size(K,1),1), K + jitter*eye(size(K)));
plot(xGrid, fSim, symbols{1}, 'linewidth',2)
xlabel('x')
ylabel('f')
title(['SE kernel - ', '$\ell$ = 1/2'], 'interpreter','latex')
set(gca,'fontsize',fontSize)

x1 = xGrid(round(nGridPoints/3));
x2 = xGrid(round(nGridPoints/2));
x3 = xGrid(round(2*nGridPoints/3));
f1 = fSim(round(nGridPoints/3));
f2 = fSim(round(nGridPoints/2));
f3 = fSim(round(2*nGridPoints/3));
ylimits = get(gca,'ylim');
xlimits = get(gca,'xlim');
ymin = ylimits(1);
ymax = ylimits(2);
xmin = xlimits(1);
xmax = xlimits(2);
plot(x1,f1,'ro')
line([x1 x1],[ymin f1],'linestyle','--','color','r')
line([xmin x1],[f1 f1],'linestyle','--','color','r')
text(x1+0.01*(xmax-xmin),f1,'$f(x_1)$','interpreter','latex')

plot(x2,f2,'ro')
line([x2 x2],[ymin f2],'linestyle','--','color','r')
line([xmin x2],[f2 f2],'linestyle','--','color','r')
text(x2+0.01*(xmax-xmin),f2,'$f(x_2)$','interpreter','latex')

plot(x3,f3,'ro')
line([x3 x3],[ymin f3],'linestyle','--','color','r')
line([xmin x3],[f3 f3],'linestyle','--','color','r')
text(x3+0.01*(xmax-xmin),f3,'$f(x_3)$','interpreter','latex')


print(strcat(figureFolder,'OneParameterAtEveryPoint'),'-dpng')


%% Smooth function, nearby points
figure('name','SmoothNearby')
subplot(1,2,1) % GP realization
KernelName = 'SquaredExp';
hyperParam = [1 1];
hyperParamNames = {'\sigma_f','\ell'}; 
plotXPoints = [20 25];
KPoints = PlotGPRealizationWithMarks(KernelName, hyperParam, hyperParamNames, plotXPoints, nGridPoints, fontSize, 'off', jitter);

subplot(1,2,2) % Bivariate normal
sigmas = sqrt(diag(KPoints));
rho = KPoints(1,2)/(sigmas(1)*sigmas(2));
PlotBivariateNormalWithMargins([0 0], sigmas, rho, 'off', 1)

print(strcat(figureFolder,'GPSmoothNearby'),'-dpng')


%% Smooth function, distant points
figure('name','SmoothDistant')
subplot(1,2,1) % GP realization
KernelName = 'SquaredExp';
hyperParam = [1 1];
hyperParamNames = {'\sigma_f','\ell'}; 
plotXPoints = [20 50];
KPoints = PlotGPRealizationWithMarks(KernelName, hyperParam, hyperParamNames, plotXPoints, nGridPoints, fontSize, 'off', jitter);

subplot(1,2,2) % Bivariate normal
sigmas = sqrt(diag(KPoints));
rho = KPoints(1,2)/(sigmas(1)*sigmas(2));
PlotBivariateNormalWithMargins([0 0], sigmas, rho, 'off', 1)

print(strcat(figureFolder,'GPSmoothDistant'),'-dpng')


%% Jagged function, nearby points
figure('name','JaggedNearby')
subplot(1,2,1) % GP realization
KernelName = 'SquaredExp';
hyperParam = [1 0.1];
hyperParamNames = {'\sigma_f','\ell'}; 
plotXPoints = [20 25];
KPoints = PlotGPRealizationWithMarks(KernelName, hyperParam, hyperParamNames, plotXPoints, nGridPoints, fontSize, 'off', jitter);

subplot(1,2,2) % Bivariate normal
sigmas = sqrt(diag(KPoints));
rho = KPoints(1,2)/(sigmas(1)*sigmas(2));
PlotBivariateNormalWithMargins([0 0], sigmas, rho, 'off', 1)

print(strcat(figureFolder,'GPJaggedNearby'),'-dpng')


%% Jagged function, distant points
figure('name','JaggedDistant')
subplot(1,2,1) % GP realization
KernelName = 'SquaredExp';
hyperParam = [1 0.1];
hyperParamNames = {'\sigma_f','\ell'}; 
plotXPoints = [20 50];
KPoints = PlotGPRealizationWithMarks(KernelName, hyperParam, hyperParamNames, plotXPoints, nGridPoints, fontSize, 'off', jitter);

subplot(1,2,2) % Bivariate normal
sigmas = sqrt(diag(KPoints));
rho = KPoints(1,2)/(sigmas(1)*sigmas(2));
PlotBivariateNormalWithMargins([0 0], sigmas, rho, 'off', 1)

print(strcat(figureFolder,'GPJaggedDistant'),'-dpng')

%% Simulating a GP

%figure('name','Simulating a GP')
%hold on
KernelName = 'SquaredExp';
hyperParam = [0.3 0.1];
hyperParamNames = {'\sigma_f','\ell'}; 
meanFunc = @(x) sin(10*x);
xGrid = linspace(0,1,1000)';
nSim = 5;
[f, K] =  GPRealizations(KernelName, hyperParam, hyperParamNames, meanFunc, xGrid, nSim, fontSize, 'on', 10^-10, 1);
set(gca,'fontsize',fontSize)
ylabel('$f(x)$')
xlabel('$x$')
print(strcat(figureFolder,'SimulateGP'),'-dpng')

%% The effect of the length scale

xGrid = linspace(-1,1,nGridPoints)';
figure('name','Realization SE and Matern')

% Square exp

K = KernelFunctions('SquaredExp', xGrid, xGrid, [1 1/10]);
subplot(2,3,1)
hold on
for i = 1:3
    fSim = mvnrnd(zeros(size(K,1),1), K + jitter*eye(size(K)));
    plot(xGrid, fSim, symbols{i}, 'linewidth',2)
end
xlabel('x')
ylabel('f')
title(['SE kernel - ', '$\ell$ = 0.1'], 'interpreter','latex')
set(gca,'fontsize',fontSize)





K = KernelFunctions('SquaredExp', xGrid, xGrid, [1 1/2]);

subplot(2,3,2)
hold on
for i = 1:3
    fSim = mvnrnd(zeros(size(K,1),1), K + jitter*eye(size(K)));
    plot(xGrid, fSim, symbols{i}, 'linewidth',2)
end
xlabel('x')
ylabel('f')
title(['SE kernel - ', '$\ell$ = 0.5'], 'interpreter','latex')
set(gca,'fontsize',fontSize)



K = KernelFunctions('SquaredExp', xGrid, xGrid, [1 1]);

subplot(2,3,3)
hold on
for i = 1:3
    fSim = mvnrnd(zeros(size(K,1),1), K + jitter*eye(size(K)));
    plot(xGrid, fSim, symbols{i}, 'linewidth',2)
end
xlabel('x')
ylabel('f')
title(['SE kernel - ', '$\ell$ = 1'], 'interpreter','latex')
set(gca,'fontsize',fontSize)




% Matern 3/2 kernel
K = KernelFunctions('Matern32', xGrid, xGrid, [1 1/10]);
subplot(2,3,4)
hold on
for i = 1:3
    fSim = mvnrnd(zeros(size(K,1),1), K + jitter*eye(size(K)));
    plot(xGrid, fSim, symbols{i}, 'linewidth',2)
end
xlabel('x')
ylabel('f')
title(['Matern32 - ', '$\ell$ = 0.1'], 'interpreter','latex')
set(gca,'fontsize',fontSize)



% Matern 3/2 kernel
K = KernelFunctions('Matern32', xGrid, xGrid, [1 1/2]);
subplot(2,3,5)
hold on
for i = 1:3
    fSim = mvnrnd(zeros(size(K,1),1), K + jitter*eye(size(K)));
    plot(xGrid, fSim, symbols{i}, 'linewidth',2)
end
xlabel('x')
ylabel('f')
title(['Matern32 - ', '$\ell$ = 0.5'], 'interpreter','latex')
set(gca,'fontsize',fontSize)




% Matern 3/2 kernel
K = KernelFunctions('Matern32', xGrid, xGrid, [1 1]);
subplot(2,3,6)
hold on
for i = 1:3
    fSim = mvnrnd(zeros(size(K,1),1), K + jitter*eye(size(K)));
    plot(xGrid, fSim, symbols{i}, 'linewidth',2)
end
xlabel('x')
ylabel('f')
title(['Matern32 - ', '$\ell$ = 1'], 'interpreter','latex')
set(gca,'fontsize',fontSize)


print(strcat(figureFolder,'RealizationsSEandMatern'),'-dpng')



%% The effect of the variance

xGrid = linspace(-1,1,nGridPoints)';
figure('name','Realization SE and Matern')

% Square exp

K = KernelFunctions('SquaredExp', xGrid, xGrid, [0.1 1/2]);
subplot(2,3,1)
hold on
for i = 1:3
    fSim = mvnrnd(zeros(size(K,1),1), K + jitter*eye(size(K)));
    plot(xGrid, fSim, symbols{i}, 'linewidth',2)
end
xlabel('x')
ylabel('f')
title(['SE kernel - ', '$\sigma_f$ = 1/10'], 'interpreter','latex')
set(gca,'fontsize',fontSize)





K = KernelFunctions('SquaredExp', xGrid, xGrid, [1 1/2]);

subplot(2,3,2)
hold on
for i = 1:3
    fSim = mvnrnd(zeros(size(K,1),1), K + jitter*eye(size(K)));
    plot(xGrid, fSim, symbols{i}, 'linewidth',2)
end
xlabel('x')
ylabel('f')
title(['SE kernel - ', '$\sigma_f$ = 1'], 'interpreter','latex')
set(gca,'fontsize',fontSize)



K = KernelFunctions('SquaredExp', xGrid, xGrid, [3 1/2]);

subplot(2,3,3)
hold on
for i = 1:3
    fSim = mvnrnd(zeros(size(K,1),1), K + jitter*eye(size(K)));
    plot(xGrid, fSim, symbols{i}, 'linewidth',2)
end
xlabel('x')
ylabel('f')
title(['SE kernel - ', '$\sigma_f$ = 3'], 'interpreter','latex')
set(gca,'fontsize',fontSize)




% Matern 3/2 kernel
K = KernelFunctions('Matern32', xGrid, xGrid, [0.1 1/2]);
subplot(2,3,4)
hold on
for i = 1:3
    fSim = mvnrnd(zeros(size(K,1),1), K + jitter*eye(size(K)));
    plot(xGrid, fSim, symbols{i}, 'linewidth',2)
end
xlabel('x')
ylabel('f')
title(['Matern32 - ', '$\sigma_f$ = 1/10'], 'interpreter','latex')
set(gca,'fontsize',fontSize)



% Matern 3/2 kernel
K = KernelFunctions('Matern32', xGrid, xGrid, [1 1/2]);
subplot(2,3,5)
hold on
for i = 1:3
    fSim = mvnrnd(zeros(size(K,1),1), K + jitter*eye(size(K)));
    plot(xGrid, fSim, symbols{i}, 'linewidth',2)
end
xlabel('x')
ylabel('f')
title(['Matern32 - ', '$\sigma_f$ = 1'], 'interpreter','latex')
set(gca,'fontsize',fontSize)




% Matern 3/2 kernel
K = KernelFunctions('Matern32', xGrid, xGrid, [3 1/2]);
subplot(2,3,6)
hold on
for i = 1:3
    fSim = mvnrnd(zeros(size(K,1),1), K + jitter*eye(size(K)));
    plot(xGrid, fSim, symbols{i}, 'linewidth',2)
end
xlabel('x')
ylabel('f')
title(['Matern32 - ', '$\sigma_f$ = 3'], 'interpreter','latex')
set(gca,'fontsize',fontSize)

% Set the same scale in all subplots
for i = 1:6
    subplot(2,3,i)
    set(gca,'ylim',[-4 4])
end

print(strcat(figureFolder,'RealizationsSEandMaternSigmaF'),'-dpng')






%% The mean can be non-zero

xGrid = linspace(-1,1,nGridPoints)';
meanFunc = sin(3*xGrid);
figure('name','Realization SE and Matern')

% Square exp

K = KernelFunctions('SquaredExp', xGrid, xGrid, [0.1 1/2]);
subplot(2,3,1)
hold on
for i = 1:3
    fSim = mvnrnd(meanFunc, K + jitter*eye(size(K)));
    plot(xGrid, fSim, symbols{i}, 'linewidth',2)
end
xlabel('x')
ylabel('f')
title(['SE kernel - ', '$\sigma_f$ = 1/10'], 'interpreter','latex')
set(gca,'fontsize',fontSize)





K = KernelFunctions('SquaredExp', xGrid, xGrid, [1/2 1/2]);

subplot(2,3,2)
hold on
for i = 1:3
    fSim = mvnrnd(meanFunc, K + jitter*eye(size(K)));
    plot(xGrid, fSim, symbols{i}, 'linewidth',2)
end
xlabel('x')
ylabel('f')
title(['SE kernel - ', '$\sigma_f$ = 1/2'], 'interpreter','latex')
set(gca,'fontsize',fontSize)



K = KernelFunctions('SquaredExp', xGrid, xGrid, [1 1/2]);

subplot(2,3,3)
hold on
for i = 1:3
    fSim = mvnrnd(meanFunc, K + jitter*eye(size(K)));
    plot(xGrid, fSim, symbols{i}, 'linewidth',2)
end
xlabel('x')
ylabel('f')
title(['SE kernel - ', '$\sigma_f$ = 1'], 'interpreter','latex')
set(gca,'fontsize',fontSize)




% Matern 3/2 kernel
K = KernelFunctions('Matern32', xGrid, xGrid, [0.1 1/2]);
subplot(2,3,4)
hold on
for i = 1:3
    fSim = mvnrnd(meanFunc, K + jitter*eye(size(K)));
    plot(xGrid, fSim, symbols{i}, 'linewidth',2)
end
xlabel('x')
ylabel('f')
title(['Matern32 - ', '$\sigma_f$ = 1/10'], 'interpreter','latex')
set(gca,'fontsize',fontSize)



% Matern 3/2 kernel
K = KernelFunctions('Matern32', xGrid, xGrid, [1/2 1/2]);
subplot(2,3,5)
hold on
for i = 1:3
    fSim = mvnrnd(meanFunc, K + jitter*eye(size(K)));
    plot(xGrid, fSim, symbols{i}, 'linewidth',2)
end
xlabel('x')
ylabel('f')
title(['Matern32 - ', '$\sigma_f$ = 1/2'], 'interpreter','latex')
set(gca,'fontsize',fontSize)




% Matern 3/2 kernel
K = KernelFunctions('Matern32', xGrid, xGrid, [1 1/2]);
subplot(2,3,6)
hold on
for i = 1:3
    fSim = mvnrnd(meanFunc, K + jitter*eye(size(K)));
    plot(xGrid, fSim, symbols{i}, 'linewidth',2)
end
xlabel('x')
ylabel('f')
title(['Matern32 - ', '$\sigma_f$ = 1'], 'interpreter','latex')
set(gca,'fontsize',fontSize)

% Set the same scale in all subplots
for i = 1:6
    subplot(2,3,i)
    set(gca,'ylim',[-3 3])
end

print(strcat(figureFolder,'RealizationsSEandMaternSineMean'),'-dpng')






%% Discrete covariates and ARD
xGrid = linspace(-1,1,nGridPoints)';
yGrid0 = 0*ones(length(xGrid),1);
yGrid1 = 1*ones(length(xGrid),1);
xGrid = [repmat(xGrid,2,1) [yGrid0;yGrid1]];

figure('name','Discrete covariates')


% Small lDiscrete
subplot(1,2,1)
K = KernelFunctions('SquaredExpARD', xGrid, xGrid, [1 1 0.1]); % Covariance function when the discrete covariate is zero
fSim = mvnrnd(zeros(size(K,1),1), K + jitter*eye(size(K)));

nGrid = size(xGrid,1);
plot3(xGrid(1:(nGrid/2),1),xGrid(1:(nGrid/2),2), fSim(1:(nGrid/2)), 'linewidth',2,'color','b')
hold on
plot3(xGrid((nGrid/2)+1:end,1),xGrid((nGrid/2)+1:end,2), fSim((nGrid/2)+1:end), 'linewidth',2,'color','r')
set(gca,'ylim',[-0.5 1.5])
set(gca,'YTick',[0 1])
zlimits = get(gca,'zlim');
plot3(xGrid(1:(nGrid/2),1),yGrid0,zlimits(1)*ones(nGrid/2,1), 'linewidth',0.5,'color', 'k', 'linestyle',':')
plot3(xGrid(1:(nGrid/2),1),yGrid1,zlimits(1)*ones(nGrid/2,1), 'linewidth',0.5,'color', 'k', 'linestyle',':')
xlabel('$x_1$', 'interpreter','latex')
ylabel('$x_2$', 'interpreter','latex')
az = -131; el = 12; view(az, el);
title('$\ell_1 = 1$, $\ell_2 = 0.1$','interpreter','latex')


% Large lDiscrete
subplot(1,2,2)
K = KernelFunctions('SquaredExpARD', xGrid, xGrid, [1 1 10]); % Covariance function when the discrete covariate is zero
fSim = mvnrnd(zeros(size(K,1),1), K + jitter*eye(size(K)));

nGrid = size(xGrid,1);
plot3(xGrid(1:(nGrid/2),1),xGrid(1:(nGrid/2),2), fSim(1:(nGrid/2)), 'linewidth',2,'color','b')
hold on
plot3(xGrid((nGrid/2)+1:end,1),xGrid((nGrid/2)+1:end,2), fSim((nGrid/2)+1:end), 'linewidth',2,'color','r')
set(gca,'ylim',[-0.5 1.5])
set(gca,'YTick',[0 1])
zlimits = get(gca,'zlim');
plot3(xGrid(1:(nGrid/2),1),yGrid0,zlimits(1)*ones(nGrid/2,1), 'linewidth',0.5,'color', 'k', 'linestyle',':')
plot3(xGrid(1:(nGrid/2),1),yGrid1,zlimits(1)*ones(nGrid/2,1), 'linewidth',0.5,'color', 'k', 'linestyle',':')
xlabel('$x_1$', 'interpreter','latex')
ylabel('$x_2$', 'interpreter','latex')
az = -131; el = 12; view(az, el);
title('$\ell_1 = 1$, $\ell_2 = 10$','interpreter','latex')


print(strcat(figureFolder,'RealizationsDiscreteCovariate'),'-dpng')

%% Squash f


% Squashing a linear function
figure('name','Squashing linear')

subplot(1,2,1)
% Linear function
x = linspace(-3,3,1000);
f = 0.2 + 0.5*x;
plot(x,f, 'linewidth',2,'color','red')
set(gca,'xlim',[-3 3])
title('Linear function')
xlabel('$x$')
ylabel('$w_0 + w_1 \cdot x$')
set(gca,'fontsize',16)

% Linear function squashed through logistic
subplot(1,2,2)
plot(x,1./(1+exp(-f)), 'linewidth',2,'color','red')
set(gca,'xlim',[-3 3])
title('Squashed linear function')
xlabel('$x$')
ylabel('$\sigma(w_0 + w_1 \cdot x)$')
set(gca,'fontsize',16)

print(strcat(figureFolder,'SquashingLinear'),'-dpng')

% Squashing a GP
figure('name','Squashing linear')

subplot(1,2,1)
% Nonlinear GP realization
KernelName = 'SquaredExp';
hyperParam = [1 1];
hyperParamNames = {'\sigma_f','\ell'}; 
meanFunc = @(x) zeros(length(x),1);
xGrid = x';
nSim = 1;
[f, K] =  GPRealizations(KernelName, hyperParam, hyperParamNames, meanFunc, xGrid, nSim, fontSize, 'on', 10^-10, 0);

plot(x,f, 'linewidth',2,'color','red')
set(gca,'xlim',[-3 3])
title('GP realization')
xlabel('$x$')
ylabel('$f(x)$')
set(gca,'fontsize',16)


% GP squashed through logistic
subplot(1,2,2)
plot(x,1./(1+exp(-f)), 'linewidth',2,'color','red')
set(gca,'xlim',[-3 3])
title('Squashed GP realization')
xlabel('$x$')
ylabel('$\sigma(f(x))$')
set(gca,'fontsize',16)

print(strcat(figureFolder,'SquashingGP'),'-dpng')