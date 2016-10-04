% Code snippet for the state space part of Advanced ML course
% Mattias Villani, http://mattiasvillani.com

%% Prelims
figureFolder = '/Users/mv/Dropbox/Teaching/AdvML/StateSpace/Slides/Graphics/';
set(0,'DefaultTextInterpreter', 'latex')
fontSize = 14;
fontSizeTitle = 16;
symbols = {'b','r','g'}; % Plot Symbols


%% Simulate AR(1) process
sigma = 1;
T = 100;
rhos = [0 0.5 0.9 1];
figure('name','AR1 realization')

for j = 1:4 
    rho = rhos(j);
    x = zeros(3*T+1,1);
    for i = 2:((3*T)+1)
        x(i) = rho*x(i-1) + sigma*randn;
    end
    x = x(2*T+2:end,1);
    subplot(2,2,j)
    if j == 1
        h = plot(1:T, x, 'o-', 'linewidth', 1);
        set(h,'markersize',2, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')
        
    else
        plot(1:T, x, 'linewidth', 2)
    end
   
    xlabel('time')
    ylabel('$y_t$')
    set(gca, 'fontsize', fontSize)
    title(['$\rho=$',num2str(rho,2)], 'fontsize',fontSizeTitle)
end
print(strcat(figureFolder,'PlotAR1'),'-dpng')


%% Simulate Two regime AR(1) process
T = 500;
sigma1 = [1 1];
sigma2 = [1 5];
rhos1 = [0.0  0.0];
rhos2 = [1.0  0.0];
TransProbs = [0.98 0.02; 0.01 0.99];
figure('name','Two-state AR1 realization')

for j = 1:2
    rho = [rhos1(j) rhos2(j)];
    sigma = [sigma1(j) sigma2(j)];
    
    % Generate the hidden state from a markov chain
    z = nan(3*T+1,1); % Start the process in state 1
    z(1) = 1;
    for i = 2:((3*T)+1)
        z(i) = 1 + (rand > TransProbs(z(i-1), 1));
    end
    
    x = zeros(3*T+1,1);
    for i = 2:((3*T)+1)
        x(i) = rho(z(i))*x(i-1) + sigma(z(i))*randn;
    end
    z = z(2*T+2:end,1);
    x = x(2*T+2:end,1);
    
    subplot(2,2,1+2*(j-1))
    plot(1:T, z, 'linewidth', 2)
    xlabel('time')
    ylabel('$x_t$')
    set(gca, 'fontsize', fontSize, 'xlim',[0 T])
    title(['$A_{11}=$',num2str(TransProbs(1,1),2), ...
        ', $A_{22}=$',num2str(TransProbs(2,2),2)], 'fontsize',fontSizeTitle)

    
    subplot(2,2,2+2*(j-1))
    plot(1:T, x, 'linewidth', 2)
    xlabel('time')
    ylabel('$y_t$')
    set(gca, 'fontsize', fontSize, 'xlim', [0 T])
    title(['$\rho_1=$',num2str(rho(1),2),', $\rho_2=$',num2str(rho(2),2), ...
        ', $\sigma_1=$',num2str(sigma(1),2),', $\sigma_2=$',num2str(sigma(2),2)],...
        'fontsize',fontSizeTitle)
end
print(strcat(figureFolder,'PlotHiddenMarkovAR1'),'-dpng')

