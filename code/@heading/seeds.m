function seeds = seeds(obj,data,vxs)
% Generate parameter seeds for the non-linear search
%
% Syntax:
%   seeds = obj.seeds(data,vxs)
%
% Description:
%   Generates a set of seed parameters for each voxel/vertex in vxs.
%
% Inputs:
%   data                  - A matrix [v t] or cell array of such
%                           matricies. The fMRI time-series data across t
%                           TRs, for v vertices / voxels. The data should
%                           have bassed through the prep stage.
%   vxs                   - Vector. A list of vertices/voxels to be
%                           processed.
%
% Optional key/value pairs:
%   none
%
% Outputs:
%   seeds                 - Cell array. Each cell contains a matrix of
%                           [v nParams] and is one of the seed sets.
%


% Derived vars
totalVxs = size(data{1},1);

% Generate default seeds
x0 = obj.initial;

% Obj variables
stimAcqGroups = obj.stimAcqGroups;
stimTime = obj.stimTime;
filterResponse = obj.filterResponse;
stimCols = size(filterResponse,2);

% Create the HRF
hrf = obj.flobsbasis*x0(end-2:end)';

% Normalize the kernel to have unit area
hrf = hrf/sum(abs(hrf));

% Create a regression matrix for the HRF implied by the FLOBS params
X = zeros(size(obj.dataTime,1),stimCols);

for ss = 1:stimCols
    
    % Grab this stimulus row
    neuralSignal = filterResponse(:,ss);
    
    % Convolve the neuralSignal by the hrf, respecting acquisition boundaries
    fit = conv2run(neuralSignal,hrf,stimAcqGroups);
    
    % If the stimTime variable is not empty, resample the fit to match
    % the temporal support of the data.
    if ~isempty(stimTime)
        dataAcqGroups = obj.dataAcqGroups;
        dataTime = obj.dataTime;
        fit = resamp2run(fit,stimAcqGroups,stimTime,dataAcqGroups,dataTime);
    end
    
    % Apply the cleaning step
    fit = obj.clean(fit);
    
    % Add the vector to the design matrix
    X(:,ss) = fit;
    
end

% Loop over vxs and generate the seed values
filterBinIdx = obj.nFixedParamsAdapt+1:obj.nFixedParamsAdapt+obj.nFilterBins;
data = catcell(2,data);
seedMatrix = repmat(x0,totalVxs,1);
for vv=1:length(vxs)

    % Obtain the time-series data and clean it
    datats = data(vv,:)';
    datats = obj.clean(datats);

    % Obtain the beta values
    b=lasso(X,datats,'Lambda',obj.lassoRegularization);

    % Store the values in the seed matrix
    seedMatrix(vxs(vv),filterBinIdx) = b';

end

% Put the seed matrix in a cell
seeds = {seedMatrix};

end


