% This code is taken with minor modifications from B. A. Surya, J. Appl. 
% Prob. 45:135-149 (2008).

% This function is analogous to Surya's LEXP function, but defined for the
% Laplace exponent required in our particular process.

% Function inputs:
% y - argument of the Laplace exponent psi(y), i.e. vector of grid points
% params - vector of model parameter values (lambda, k, theta)

% Function output:
% psi - the Laplace exponent corresponding to a gamma-distributed lifetime
% in the age-dependent extinction model


function psi = LEXP(y,params)

% extract required model parameters:
lam = params(1); % speciation rate
k = params(2); % shape param of the lifetime dist
theta = params(3); % scale param of the lifetime dist

% Evaluating the closed-form expression for psi in the case of a gamma
% lifetime distribution:

integr = 1 - (1+y*theta).^(-k);
psi = y - lam*integr;

end