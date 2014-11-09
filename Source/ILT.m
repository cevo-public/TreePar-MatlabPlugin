% This code is taken with minor modifications from B. A. Surya, J. Appl. 
% Prob. 45:135-149 (2008).

% This function takes the numerical inverse Laplace transform (ILT) of a  
% specified function. Surya denotes the function G=ILT(F,X,P) where G is 
% the result of numerically inverting the univariate Laplace transform 
% given by F, and P denotes any parameters of F. I have left out P in the 
% function inputs, as these parameters are instead declared globally.

% Function inputs:
% F - a function defining the Laplace transform to be inverted.
% X - the points at which to evaluate the ILT. 

% Function output:
% G - the inverse Laplace transform of F.


function G = ILT(F,X)

% Set parameters of the ILT:
    % I have left these at the same values used by Surya and obtained
    % acceptable results. However these values can be modified to control
    % error: the absolute error bound (i.e. diff btwn true and numerical 
    % ILT) is given by Surya's equation (17): C*exp(-A1)/(1-exp(-A1)), 
    % where in our notation (as in manuscript), 
    % C = 1/psi'(eta) = 1 - lambda*int_0^Inf(pi(dx) x exp(-eta*x)), 
    % and eta is the largest root of psi.

N=11; M=9; A1 = 14.0; l1 = 1;

% Create weights to be used in the Euler summation of the partial sums:

mx = pascal(M+1);
my = fliplr(mx);
bn = diag(my)*2^(-M);
weight = ones([2*N+1 1]);
head = cumsum(bn);
tail = 1-cumsum(bn);
tail(M+1,:) = []; % don't understand why this is needed - makes it a mtx??
head(M+1,:) = [];
weight = [head; weight; tail];

% Set values of arguments at which transform series is to be evaluated:

val1 = -(N+M):(N+M);
val1 = (1i*pi*val1 + A1/2)/l1;
X_inv = 1./X; % note: not defined (infinite) @ x=0
X_args = kron(X_inv,val1); %?

% Evaulate the integrand at all the points:

integrand = feval(F,X_args).*exp(X_args*diag(kron(X,ones(1,1+2*N+2*M))));

% Prepare the matrix which will post-multiply the integrand:

right = kron(diag(X_inv),weight)/(2*l1);

% Finally, the inverse Laplace transform of F is given by:

G = real(integrand*right);

end