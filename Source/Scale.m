% This code is taken with minor modifications from B. A. Surya, J. Appl. 
% Prob. 45:135-149 (2008).

% This function computes the scale function W(x) required in the 
% computation of the likelihood, by combining Surya's functions 'Scale' and
%'qScale' for the special case q=0 which we require.

% Function inputs:
% x - argument of the scale function W(x), i.e. vector of grid points on 
% which the scale function is to be evaluated. x values should be strictly
% greater than zero for the ILT to be defined!
% params - vector of model parameter values (lambda, k, theta)

% Function output:
% W - the scale function W(x) required for the likelihood computation, at 
% the given model parameter values.


function W = Scale(x,params)

% First compute the inverse Laplace transform of the function 'funcscale'
% evaluated at the points x. Wtemp corresponds to W_Phi(q) (x) in Surya's 
% notation, where Phi(q) = eta for us.

Wtemp = ILT(@(y) funcscale(y,params),x);

% Then compute the actual scale function we want, W, which corresponds to 
% W^(q) (x) with q=0 in Surya's notation.  This makes use of Surya's 
% Proposition 1.

W = exp(LEXProot(params)*x).*Wtemp;

end