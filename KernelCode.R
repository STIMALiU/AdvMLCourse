# Matern32  kernel
k <- function(sigmaf = 1, ell = 1)  
{   
	rval <- function(x, y = NULL) 
	{	r = sqrt(crossprod(x-y))
		 return(sigmaf^2*(1+sqrt(3)*r/ell)*exp(-sqrt(3)*r/ell))   
	}   
	class(rval) <- "kernel"   
	return(rval) 
} 