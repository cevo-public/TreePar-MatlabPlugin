% Function to maximize likelihood starting from a given initial point.
% Outputs maximum likelihood estimates (parameter values) and negative log
% likelihood at this point. Requires that the data (tree branching times) 
% has already been read in by ReadTreeFcn and loads the mat-file output 
% from ReadTreeFcn.



% INPUTS:

% setupFilename: name of mat file (w/o extension) containing intermediate 
% variables (produced by ReadTreeFcn)

% outputFilename: name of plain text file (w/o extension) where results of 
% inference (MLEs and maximized log likelihood value) will be written

% outputPrec: number of decimal places to print output

% Lcond: 'C' for conditioning likelihood on crown age, 'S' for conditioning
% on stem age

% model: 'E' to run inference under exponential model, 'G' under gamma model

% Init: initial point from which to run optimization, in form 'lam k
% th' for gamma model or 'lam th' for expmodel

% p: sampling probability



% Note that this will be run from a shell script, so all inputs are taken
% as strings. str2num can convert vectors and will interpret e.g. '[lam k th]' 
% as vector of length 3.



function MaxLFcn(setupFilename,outputFilename,outputPrec,Lcond,model,Init,p)


outputPrec = str2num(outputPrec);
Init = str2num(Init);
p = str2num(p);


% Load tree-related variables set by ReadTree:
load(setupFilename);



%%%%%%%%%% RUN MAX LIKELIHOOD INFERENCE %%%%%%%%%%

% Select optimization algorithm:
options = optimset('Algorithm','interior-point');


if(strcmpi(model,'E'))
    % runs max llhd inference under exponentially-distributed lifetime


    % Set up constraints to use in optimization (see fmincon documentation):
    % here we constrain that our parameters (lambda, theta) are non-negative
    A = diag(-ones(1,2));
    b = zeros(2,1);


    % Call inference function, inputting variables defined above. 
    % Use likelihood conditioned on crown age or on stem age, depending on 
    % how Lcond is set.
    
    if(strcmpi(Lcond,'C')) % conditioning on crown age of tree
        
        [MLEpars,fval] = fmincon(@(z) LikeliCrown(x,[z(1) 1 z(2) p], Ctpts), ...
                        Init', A, b, [],[],[],[],[],options);
                 
    elseif(strcmpi(Lcond,'S')) % conditioning on stem age of tree
        
        [MLEpars,fval] = fmincon(@(z) LikeliStem(x,[z(1) 1 z(2) p], Ctpts), ...
                        Init', A, b, [],[],[],[],[],options);
        
    end



elseif(strcmpi(model,'G'))
    % runs max llhd inference under gamma-distributed lifetime
    
    
    % Set up constraints to use in optimization (see fmincon documentation):
    % here we constrain that our parameters (lambda, k, theta) are non-negative
    A = diag(-ones(1,3));
    b = zeros(3,1);

    
    % Call inference function, inputting variables defined above. 
    % Use likelihood conditioned on crown age or on stem age, depending on 
    % how Lcond is set.
    
    if(strcmpi(Lcond,'C')) % conditioning on crown age of tree
        
        [MLEpars fval] = fmincon(@(z) LikeliCrown(x,[z(1) z(2) z(3) p], Ctpts),...
                           Init', A, b, [],[],[],[],[],options);
        
    elseif(strcmpi(Lcond,'S')) % conditioning on stem age of tree
        
        [MLEpars fval] = fmincon(@(z) LikeliStem(x,[z(1) z(2) z(3) p], Ctpts),...
                           Init', A, b, [],[],[],[],[],options);
    end
    
    
end  % end if/elseif inference




%%%%%%%%%% OUTPUT %%%%%%%%%%

precString = ['%.' num2str(outputPrec) 'f'];
dlmwrite([outputFilename '.txt'],[MLEpars' fval],'precision',precString);



end