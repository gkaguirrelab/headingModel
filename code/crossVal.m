function crossVal(sub, roi_mask)
% crossValDemo.m
% tbUseProject('headingModel')
% This script demonstrates working with model objects for the purpose of
% performing cross-validation analyses.


% The mat file we need is located in the "data" directory. Construct a path
% to this file, starting with the location of this script.
% parcels = {'bMask', 'EVC', 'OPA', 'PPA', 'RSC', 'PHG', 'ERC', 'Hipp', 'Thal'...
%     'VTC'};
% sub = 'sub-08';
% roi_mask='RSC';
train_dst = 'A';
test_dst = 'B';
nFitBins = 45;
n_par_workers=4;
% The TR of the experiment, in seconds
tr=2;
% Define modelOpts
hrfSearch = true;
adaptSearch = true;
lassoRegularization=0.05; %1; 0.05; 0.25; 0.5;
penalty_degree='005';

modelOpts = {'nFilterBins',nFitBins,'hrfSearch',hrfSearch,'adaptSearch',adaptSearch,'typicalGain',1, ...
'lassoRegularization',lassoRegularization};
parpool(n_par_workers);

% Load the fMRI data
fileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),['data/' sub], ...
[sub '_city1' train_dst '_stimulus_data_' roi_mask '.mat']);
load(fileName,'data', 'stimulus')
data_train=data;
stimulus_train_og=stimulus;

nDTs_og = size(stimulus_train_og{1},2);
nTRs = size(data_train{1},2);
deltaTStim_og = 2*nTRs/nDTs_og; % The temporal sampling interval (in seconds) of the
                       % heading direction vector
% Resample the heading direction vector
headingDownsampleFactor = 3;
stimulus_train = cellfun(@(x) wrapTo2Pi(decimate(unwrap(x),headingDownsampleFactor)),stimulus_train_og,'UniformOutput',false);
% nDTs = length(stimulus_train{1});
deltaTStim = deltaTStim_og * headingDownsampleFactor;

% How many samples do we need to pad the stimulus to have the requested
% duration of padding in seconds?
nTRsToPad = 8;
nDTsToPad = round(2*nTRsToPad*(1/deltaTStim));
time_pad='8TR-pad';

% Extend the stimulus back in time to provide some "warm up" for the HRF;
% create a stimTime array
train_stimTime={};
for ii=1:length(stimulus_train)

    % Grab an acquisition from the cell array
    thisStim = stimulus_train{ii};

    % Create a stimTime vector
    train_stimTime{ii} = -nDTsToPad*deltaTStim:deltaTStim:deltaTStim*(length(thisStim)-1);

    % Pad this stim
    thisStim = [repmat(thisStim(1),1,nDTsToPad), thisStim];
    stimulus_train{ii} = thisStim;

end

% Define the set of vertices
vxs = 1:size(data_train{1},1);
% vxs = 1;
% Call the forwardModel
% results = forwardModel(data,stimulus,tr,'modelClass','heading','vxs',vxs,'modelOpts',modelOpts);
results = forwardModel(data_train,stimulus_train,tr,'stimTime',train_stimTime,'modelClass','heading','vxs',vxs,'modelOpts',modelOpts);

%% VALIDATION
% Note that the steps below will over-write the "stimulus" and "data"
% variables. We may wish to make this more robust later.

% Load the validation dataset
fileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),['data/' sub],...
    [sub '_city1' test_dst '_stimulus_data_' roi_mask '.mat']);
load(fileName,'data', 'stimulus')
data_test = data;
stimulus_test_og = stimulus;
stimulus_test = cellfun(@(x) wrapTo2Pi(decimate(unwrap(x),headingDownsampleFactor)),stimulus_test_og,'UniformOutput',false);
% Extend the stimulus back in time to provide some "warm up" for the HRF;
% create a stimTime array
test_stimTime={};
for ii=1:length(stimulus_test)

    % Grab an acquisition from the cell array
    thisStim = stimulus_test{ii};

    % Create a stimTime vector
    test_stimTime{ii} = -nDTsToPad*deltaTStim:deltaTStim:deltaTStim*(length(thisStim)-1);

    % Pad this stim
    thisStim = [repmat(thisStim(1),1,nDTsToPad), thisStim];
    stimulus_test{ii} = thisStim;

end
% Re-create the model, now with the new stimulus and data
% model = heading(data,stimulus,tr,modelOpts{:});
model = heading(data_test,stimulus_test,tr,'stimTime',test_stimTime,modelOpts{:});

% Prepare the validation data into a signal. You would want to loop over
% vertices to test this across the entire dataset.
dataPrep = catcell(2,model.prep(data_test));
cv_R2=zeros(length(vxs), 1);
cv_R=zeros(length(vxs), 1);
parfor ii=1:length(vxs)
    signal = model.clean(dataPrep(vxs(ii),:)');
    % Obtain the parameter values of the previous fit for this vertex
    x = results.params(vxs(ii),:);
    % Obtain the R2 metric for the signal, given the params (x) we have
    % previously calculated
    valR2 = model.metric(signal,x);
    valR = calccorrelation(signal, model.forward(x));
    cv_R2(ii, 1) = valR2;
    cv_R(ii, 1) = valR;
end

% save R2 and model params
model_params=results.params;
outFileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),'results',...
    [sub '_city1' train_dst '_head-' num2str(nFitBins) ...
    '_city1' test_dst '_' roi_mask '_nl-lasso-' penalty_degree '_' time_pad '.mat']);
disp(outFileName);
save(outFileName, 'cv_R', 'cv_R2', 'model_params');

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolsize = 0;
else
    poolsize = poolobj.NumWorkers;
end
disp([num2str(poolsize) ' workers are running'])

end
