function K = PlotGPRealizationWithMarks(KernelName, hyperParam, hyperParamNames, plotXPoints, nGridPoints, fontSize, ticks, jitter)

%% Plots a GP realization with marked out points

xGrid = linspace(-1,1,nGridPoints)';

% Simulate GP
K = KernelFunctions(KernelName, xGrid, xGrid, hyperParam);
fSim = mvnrnd(zeros(size(K,1),1), K + jitter*eye(size(K)));

plot(xGrid, fSim, 'b', 'linewidth',2)
hold on
xlabel('')
ylabel('')

ylimits = get(gca,'ylim');
xlimits = get(gca,'xlim');
ymin = ylimits(1);
ymax = ylimits(2);
xmin = xlimits(1);
xmax = xlimits(2);

TitleString = '';
for i = 1:length(hyperParam)
    TitleString = strcat(TitleString,'$',hyperParamNames{i},' = ',num2str(hyperParam(i),3),'$, ');
end
title(TitleString(1:(end-1)), 'interpreter','latex')
set(gca,'fontsize',fontSize)



% Adding text to function values
for j = 1:length(plotXPoints)
    x = xGrid(plotXPoints(j));
    f = fSim(plotXPoints(j));
    plot(x,f,'ro')
    line([x x],[ymin f],'linestyle','--','color','r')
    line([xmin x],[f f],'linestyle','--','color','r')
    text(xmin-0.15*(xmax-xmin),f-0.00*(ymax-ymin),['$f(x_',int2str(j),')$'],'interpreter','latex')
    text(x-0.02*(xmax-xmin),ymin-0.02*(ymax-ymin),['$x_',int2str(j),'$'],'interpreter','latex')
end

if strcmpi(ticks,'off')
    set(gca,'YTick',[])
    set(gca,'XTick',[])
end

if isempty(plotXPoints)
    K = [];
else
    K = K(plotXPoints,plotXPoints);
end