% This code is taken with minor modifications from B. A. Surya, J. Appl. 
% Prob. 45:135-149 (2008).

% This function defines the Laplace transform function to be inverted.
% Here I have used Surya's 'funcscale' specialized to the case q=0.

% Function inputs:
% y - argument of the Laplace transform function, i.e. vector of grid points
% params - vector of model parameter values (lambda, k, theta)

% Function output:
% LTfunc - the Laplace transform function to be inverted, as a function of
% y


function LTfunc = funcscale(y,params)

% Define the Laplace transform function to invert, namely 1/psi(y+eta) for
% our purposes.
LTfunc = 1./LEXP(y+LEXProot(params),params);

end