% crossValDemo.m
% tbUseProject('headingModel')
% This script demonstrates working with model objects for the purpose of
% performing cross-validation analyses.


% The mat file we need is located in the "data" directory. Construct a path
% to this file, starting with the location of this script.
train_dst = 'A';
test_dst = 'B';
nFitBins = 45;
n_par_workers=4;
% The TR of the experiment, in seconds
tr=2;
% Define modelOpts
hrfSearch = true;
adaptSearch = false;
lassoRegularization=0;
penalty_degree='02';
nTRsToPad = 8;
time_pad='8TR-pad';
if ~adaptSearch
    time_pad='8TR-pad_adapt-0';
end
modelOpts = {'nFilterBins',nFitBins,'hrfSearch',hrfSearch,'adaptSearch',adaptSearch,'typicalGain',1, ...
'lassoRegularization',lassoRegularization};
if isempty(gcp('nocreate')) && n_par_workers>1; parpool(n_par_workers); end
% parpool(4);

sub='sub-08';
parcel='bMask';
fileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),['data/' sub], ...
[sub '_city1' train_dst '_stimulus_data_' parcel '.mat']);

% Load the stimulus and data variables
load(fileName,'stimulus','data')

% Extend the stimulus back in time to provide some "warm up" for the HRF;
% create a stimTime array
for ii=1:length(stimulus)

    % Grab an acquisition from the cell array
    thisStim = stimulus{ii};

    % Create a stimTime vector
    train_stimTime{ii} = -nTRsToPad*tr:tr:tr*(length(thisStim)-1);

    % Pad this stim
    thisStim = [repmat(thisStim(1),1,nTRsToPad), thisStim];
    stimulus{ii} = thisStim;
end

% Define the set of vertices
vxs = 1:size(data{1},1);

% Let's just work with one vertex for the demo
% ii = 48044; ii = 114371;
% vxs = vxs(ii);

% Call the forwardModel
results = forwardModel(data,stimulus,tr,'stimTime',train_stimTime,'modelClass','heading','vxs',vxs,'modelOpts',modelOpts);

%% VALIDATION
% Note that the steps below will over-write the "stimulus" and "data"
% variables. We may wish to make this more robust later.

% Load the validation dataset
fileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),['data/' sub],...
    [sub '_city1' test_dst '_stimulus_data_' parcel '.mat']);
load(fileName,'stimulus','data')
% Extend the stimulus back in time to provide some "warm up" for the HRF;
% create a stimTime array
for ii=1:length(stimulus)

    % Grab an acquisition from the cell array
    thisStim = stimulus{ii};

    % Create a stimTime vector
    test_stimTime{ii} = -nTRsToPad*tr:tr:tr*(length(thisStim)-1);

    % Pad this stim
    thisStim = [repmat(thisStim(1),1,nTRsToPad), thisStim];
    stimulus{ii} = thisStim;

end
% Re-create the model, now with the new stimulus and data
model = heading(data,stimulus,tr,'stimTime',test_stimTime,modelOpts{:});
% Prepare the validation data into a signal. You would want to loop over
% vertices to test this across the entire dataset.
dataPrep = catcell(2,model.prep(data));
% ii = 114371; % Just do the first vertex
cv_R2=zeros(length(vxs), 1);

parfor ii=1:length(vxs)
%         for ii=1:length(vxs)
    signal = model.clean(dataPrep(vxs(ii),:)');

    % Obtain the parameter values of the previous fit for this vertex
    x = results.params(vxs(ii),:);

    % Obtain the R2 metric for the signal, given the params (x) we have
    % previously calculated
    valR2 = model.metric(signal,x);
    cv_R2(ii, 1)=valR2;
end

% save R2
model_params=results.params;
outFileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),'results',...
    [sub '_city1' train_dst '_head-' num2str(nFitBins) ...
    '_city1' test_dst '_' parcel  '_nl-lasso-' penalty_degree '_' time_pad '.mat']);
save(outFileName,'cv_R2', 'model_params','results');
