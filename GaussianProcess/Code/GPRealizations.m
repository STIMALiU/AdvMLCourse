function [fSim, K] = GPRealizations(KernelName, hyperParam, hyperParamNames, meanFunc, xGrid, nSim, fontSize, ticks, jitter, doPlot)

%% Plots a GP realization with marked out points

plotSymbols = {'b','r','g','m','c'};
cmp = colormap(parula(nSim));
% Simulate GP

mu = meanFunc(xGrid);
K = KernelFunctions(KernelName, xGrid, xGrid, hyperParam);
fSim = nan(size(xGrid,1),nSim);
if doPlot
    figure('name','GP realizations')
    hold on
    stdev = sqrt(diag(K));
    lowerBand = mu - 1.96*stdev;
    upperBand = mu + 1.96*stdev;
    
    colorband = 0.9*[1 1 1];
    patchHandle = patch([xGrid' fliplr(xGrid')],[lowerBand' fliplr(upperBand')], colorband);
    set(patchHandle,'faceLighting','phong','facealpha',0.5,'edgecolor',min([1.05*colorband],[1 1 1]),'edgealpha',0.0)
    plot(xGrid, mu,'k--')
    legend({'95% Probability bands','GP mean'})
end


for i = 1:nSim
    fSim(:,i) = mvnrnd(mu, K + jitter*eye(size(K)));
    if doPlot
        plot(xGrid, fSim(:,i), 'linewidth',2)
    end
end

if doPlot
    TitleString = '';
    for i = 1:length(hyperParam)
        TitleString = strcat(TitleString,'$',hyperParamNames{i},' = ',num2str(hyperParam(i),3),'$, ');
    end
    title(TitleString(1:(end-1)), 'interpreter','latex')
    set(gca,'fontsize',fontSize)
    
    if strcmpi(ticks,'off')
        set(gca,'YTick',[])
        set(gca,'XTick',[])
    end
end
