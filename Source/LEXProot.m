% This code is taken with minor modifications from B. A. Surya, J. Appl. 
% Prob. 45:135-149 (2008).

% This function finds the largest root, eta, of the Laplace exponent, 
% psi(y), defined in the function LEXP.  This is essentially Surya's Phir 
% function specialized to the case q=0 which we require.

% Function inputs:
% params - vector of model parameter values (lambda, k, theta)
% bounds (OPTIONAL INPUT) - 2-element vector giving range in which to
% search for a root. If not specified, the search for a root begins from
% the initial point y=1 and takes a maximum of 100 steps.

% In the coding scheme here I only ever call LEXProot without bounds. 
% However, if root-finding appears problematic for your requirements, there 
% is an opportunity for troubleshooting by calling LEXProot with bounds. 
% There is also an opportunity to try to reduce numerical error by 
% increasing 'MaxFunEvals' beyond 100 if the root does not seem to be
% sufficiently accurate.

% Function output:
% eta - the largest root of the Laplace exponent psi, also representing the
% Malthusian parameter in our model.


function eta = LEXProot(params,bounds)

if nargin<2
    eta = fsolve(@(y) LEXP(y,params),1,optimset('MaxFunEvals',100,'Display','off'));
else
    eta = fzero(@(y) LEXP(y,params), bounds,optimset('Display','off'));
end

end