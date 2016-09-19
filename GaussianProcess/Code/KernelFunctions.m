function K = KernelFunctions(kernelFuncName, x, xStar, hyperParam)
% Kernel functions


switch lower(kernelFuncName)
    case 'squaredexp'
        kernelFunc = @(x,xs) hyperParam(1)^2*exp(-0.5*(abs(x-xs)^2/hyperParam(2)^2));
    case 'gaussian'
    case 'matern'
        kernelFunc = @(x,xs) hyperParam(1)^2 *(  (2^(1-hyperParam(3))/gamma(hyperParam(3))) * ...
            (sqrt(2*hyperParam(3))*abs(x-xs)/hyperParam(2))^hyperParam(3) * ...
            besselk(hyperParam(3), sqrt(2*hyperParam(3))*abs(x-xs)/hyperParam(2))   ); 
    case 'matern32'
        kernelFunc = @(x,xs) hyperParam(1)^2*(...
            (1 + sqrt(3)* abs(x-xs)/hyperParam(2) ) *...
            exp(- sqrt(3)*abs(x-xs)/hyperParam(2) )    );
    case 'matern52'
        kernelFunc = @(x,xs) hyperParam(1)^2*(...
            (1 + sqrt(5)* abs(x-xs)/hyperParam(2) + 5*abs(x-xs)^2/(2*hyperParam(2)^2) ) *...
            exp(- sqrt(5)*abs(x-xs)/hyperParam(2) )    );
    case 'squaredexpard'
        H = diag(1./hyperParam(2:end));
        %Quad = (x-xs)*H*(x-xs)';
        kernelFunc = @(x,xs) hyperParam(1)^2*exp( -0.5*( (x-xs)'*H*(x-xs) ) );

end

K = nan(length(x),length(xStar));
for i = 1:length(x)
    for j = 1:length(xStar)
        K(i,j) = kernelFunc(x(i,:)', xStar(j,:)');
    end
end
