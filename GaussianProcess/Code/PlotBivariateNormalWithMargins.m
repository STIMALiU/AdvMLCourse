
function []  = PlotBivariateNormalWithMargins(mu, sigmas, rho, ticks, titleFlag)
%% Script that plots contours of the Gaussian density in two dims

nGrid = 100;
Sigma = [sigmas(1)^2 rho*prod(sigmas);rho*prod(sigmas) sigmas(2)^1];
       
x1 = mu(1) + linspace(-3,3,nGrid);
x2 = mu(2) + linspace(-3,3,nGrid);

[X1Grid, X2Grid] = meshgrid(x1,x2);
X = [X1Grid(:) X2Grid(:)];
y = mvnpdf(X, mu, Sigma);
Y = reshape(y,nGrid,nGrid);

if titleFlag
    title(['Correlation coefficient = ',num2str(rho,2)],'fontsize',12)
end

xlabel('')
ylabel('')
hold on
posVector = get(gca,'position');
a = gca;

if strcmpi(ticks,'off')
    set(a,'YTick',[])
    set(a,'XTick',[])
end

% Plotting the marginal distribution of x1
a1 = axes;
hold on
set(a1,'position',posVector,'box','off','color','none','xtick',[],'ytick',[],'ylim',[0 1])
margDensity = 0.005+normpdf(x1,mu(1),sqrt(Sigma(1,1)));
colorband = 0.7*[1 1 1];
patchHandle = patch([x1' fliplr(x1')],[margDensity' fliplr(margDensity')],colorband);
set(patchHandle,'faceLighting','phong','facealpha',0.5,'edgecolor',min([1.05*colorband],[1 1 1]),'edgealpha',0.0)

% Plotting the marginal distribution of x2
a2 = axes;
hold on
set(a2,'position',posVector,'box','off','color','none','xtick',[],'ytick',[],'xlim',[0.02 1],'ylim',[-4 2])
margDensity = normpdf(x2,mu(2),sqrt(Sigma(2,2)));
colorband = 0.7*[1 1 1];
patchHandle = patch([margDensity' fliplr(margDensity')],[x2' fliplr(x2')],colorband);
set(patchHandle,'faceLighting','phong','facealpha',0.5,'edgecolor',min([1.05*colorband],[1 1 1]),'edgealpha',0.0)

a3 = axes;
hold on
set(a3,'position',posVector,'box','off','color','none','xtick',[],'ytick',[])
contour(X1Grid,X2Grid,Y,5,'k');

box on

set(gcf,'renderer','painters')
