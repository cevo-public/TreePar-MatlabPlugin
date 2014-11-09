% Function to evaluate the negative log likelihood at a given point (parameter
% values). Requires that the data (tree branching times) has already been
% read in by ReadTreeFcn and loads the mat-file output from ReadTreeFcn.


% INPUTS:

% setupFilename: name of mat file (w/o extension) containing intermediate 
% variables (produced by ReadTreeFcn)

% outputFilename: name of plain text file (w/o extension) where the 
% computed log likelihood value will be written

% outputPrec: number of decimal places to print output

% Lcond: 'C' for conditioning likelihood on crown age, 'S' for conditioning
% on stem age

% Param = (lam, k, th): point at which to evaluate likelihood

% p: sampling probability



% Note that this will be run from a shell script, so all inputs are taken
% as strings



function EvalLFcn(setupFilename,outputFilename,outputPrec,Lcond,Param,p)


outputPrec = str2num(outputPrec);
Param = str2num(Param);
p = str2num(p);


% Load tree-related variables:
load([setupFilename '.mat'])


%%%%%%%% CALCULATE LOG LIKELIHOOD %%%%%%%%
    
if(strcmpi(Lcond,'C'))
    neglogL = LikeliCrown(x, [Param p], Ctpts);
elseif(strcmpi(Lcond,'S'))
    neglogL = LikeliStem(x, [Param p], Ctpts);
end



%%%%%%%%%% OUTPUT %%%%%%%%%%

precString = ['%.' num2str(outputPrec) 'f'];
dlmwrite([outputFilename '.txt'],neglogL,'precision',precString)


end