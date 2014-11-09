% This function computes the negative log likelihood of a given tree (i.e. 
% set of branching times) for given model parameters, conditioned on 
% sampling at least one tip and on the CROWN AGE of the tree.

% Function inputs:
% x - vector of grid points on which to evaluate the scale function W(x)
% (assumed to be an evenly spaced grid ranging from ZERO TO THE CROWN AGE  
% of the tree)
% params - vector of model parameter values (lambda, k, theta, p)
% cTimePts - vector containing indices of x at which branching times occur

% Function output:
% neglogLlhd - negative log likelihood of the data (branching times) given
% the model parameter values

function neglogLlhd = LikeliCrown(x,params,cTimePts)

% extract sampling probability (p) from the list of model parameters
p = params(4);

% determine the grid spacing, assuming evenly spaced
diffx = x(2)-x(1);

% numerically compute the scale function W(x) given model params
W = Scale(x(2:end),params); % removes x=0 where methods are ill-defined
W = [1 W]; % adds back in the point W(0) which is known to =1

% numerically compute the derivative of W w.r.t. x: error increases with
% grid spacing diffx   
diffW2 = W(3:end) - W(1:(end-2));
Wderiv = diffW2/(2*diffx); % 2nd-order approximation for derivative on interior points
Wderiv = [(W(2)-W(1))/diffx Wderiv (W(end)-W(end-1))/diffx]; % 1st-order approx at endpoints

% COMPUTE NEGATIVE LOG LIKELIHOOD:
neglogLlhd = 2*log(1-p+p*W(end)) - sum(log(p*Wderiv(cTimePts))) ...
    + 2*sum(log(1-p+p*W(cTimePts)));

end