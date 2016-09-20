% 1.Generate example data.
rng(0,'twister');
N = 1000;
X = linspace(0,1,N)';
X = [X,X.^2];
y = 1 + X*[1;2] + sin(20*X*[1;-2]) + 0.2*randn(N,1);

% 2. Fit the model using 'SR' and predict using 'FIC'. Use 20
% points in the active set selected using 'Entropy'.
gpr = fitrgp(X,y,'KernelFunction','SquaredExponential','FitMethod','SR','PredictMethod','FIC',...
    'Basis','None','Optimizer','fminsearch','KernelParameters',[1;1],'ActiveSetSize',20,'ActiveSetMethod','Entropy',...
    'Sigma',1,'Standardize',true,'verbose',1);

% 3. Plot fit and prediction intervals.
[pred,se,ci] = predict(gpr,X,'Alpha',0.01);
figure;
plot(y,'r');
hold on;
plot(pred,'b')
plot(ci(:,1),'g--');
plot(ci(:,2),'k--');
legend('Data','Pred','Lower 95%','Upper 95%','Location','Best');



