using Distributions
using Gadfly

function SquaredExpKernel(x,xStar,param)
	# Squared exponential kernel
	sigma,lenScale = param
	sigma*exp(-0.5*((x-xStar)/lenScale)^2)
end

function Kernel2Cov(x,xStar, kernelFunc, param)
	# Function that computes n-by-n covariance matrix from kernelFunc with hyperparameters given by param
	K = nans(length(x),length(xStar))
	for (i,xi) in enumerate(x)
		for (j,xStari) in enumerate(xStar)
			K[i,j] = kernelFunc(xi,xStari,param)	
		end
	end
    return(K)
end

function MeanFunc(meanType)
	if meanType == "zero"
		return(zeros)
	end
	if meanType == "sin"
		return(sin)
	end
	if meanType == "cos"
		return(cos)
	end
end

function SimGP(x, meanType, kernelFunc, kernelParam, nSim)
  # Simulates nSim realizations (function) form a Gaussian process with mean m(x) and covariance K(x,x')
  # over a grid of inputs (x)
  n = length(x)
  f = zeros(n,nSim)
  meanVector = MeanFunc(meanType)(x)
  covMat = Kernel2Cov(x,x, kernelFunc, kernelParam)
  f = rand(MvNormal(meanVector, covMat),n) # the nSim columns of f store the individual vector draws
  return(f)
end

function PostGPRegr(y, X, Xstar, meanType, kernelFunc, kernelParam, sigmaNoise)
	# Computes the posterior distribution of f | y, X, Xstar
	n = size(X)[1]
	mx = MeanFunc(meanType)(X)
	ms = MeanFunc(meanType)(Xstar)
	Kss = Kernel2Cov(Xstar, Xstar, kernelFunc, kernelParam)
	Kxs = Kernel2Cov(X, Xstar, kernelFunc, kernelParam)
	Ksx = Kernel2Cov(Xstar, X, kernelFunc, kernelParam)
	Kxx = Kernel2Cov(X, X, kernelFunc, kernelParam)
	postCov = Kss - Ksx*(Kxx+sigmaNoise*eye(n))\Kxs
	postMean = ms + Ksx*(Kxx+sigmaNoise*eye(n))\(y-mx)
	return(postMean,postCov)
end



