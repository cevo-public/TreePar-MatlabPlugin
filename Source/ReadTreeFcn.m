% Function to read in a list of tree branching times, define intermediate 
% variables necessary to evaluate likelihood, and save these in a mat file.


% INPUTS:

% dataFilename: name of plain text file (w/o extension) containing a list 
% of branching times in the tree. These times do not have to be in any 
% particular order, but the *maximum* value will be taken in subsequent 
% likelihood evaluation to represent either the stem age or the crown age 
% of the tree, depending on user's choice of conditioning.

% setupFilename: name of mat file (w/o extension) where 
% intermediate variables will be stored

% numgridpts: number of grid points (zero to tree age) for numerical
% evaluation. The more points, the smaller the numerical error but the 
% slower the computation. Note that an *evenly spaced* grid of time points 
% from zero to (crown/stem) age of tree will be created.



% Note that this will be run from a shell script, so all inputs are taken
% as strings



function ReadTreeFcn(dataFilename,setupFilename,numgridpts)


numgridpts = str2num(numgridpts);


% Define tree age and branching times from data file:

times = dlmread([dataFilename '.txt']); % read in branching times in data file specified by user
[treeAge ind] = max(times); % treeAge will be interpreted as crown age if isCrown=1; otherwise as stem age
ctimes = times([1:(ind-1) (ind+1):end]); % the remaining branching times


% Define grid of timepoints on which to evaluate scale function:

diffx = treeAge/(numgridpts-1); % grid spacing, assumed even
x = 0:diffx:treeAge; % grid of timepoints


% Find indices of timepoint grid at which branching events occur

Ctpts = zeros(1,length(ctimes)); % initialize 

for i=1:length(ctimes) 
    [minval,tpt] = min(abs(ctimes(i)-x)); % find grid point closest to branching time
    Ctpts(i) = tpt; 
end


% save intermediate variables for use in likelihood evaluation
save([setupFilename '.mat'])


end