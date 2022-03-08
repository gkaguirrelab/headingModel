% crossValDemo.m
%
% This script demonstrates working with model objects for the purpose of
% performing cross-validation analyses.


% The mat file we need is located in the "data" directory. Construct a path
% to this file, starting with the location of this script.
sub='sub-08';
fileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),'data',[sub '_city1A_stimulus_data_bMask-city1AB-100.mat']);

% Load the stimulus and data variables
load(fileName,'stimulus','data')

% The TR of the experiment, in seconds
tr=2;

% Define the set of vertices
vxs = 1:size(data{1},1);

% Let's just work with one vertex for the demo
ii = 1;
vxs = vxs(ii);

% Define modelOpts
modelOpts = {'nFilterBins',8};

% Call the forwardModel
results = forwardModel(data,stimulus,tr,'modelClass','heading','vxs',vxs,'modelOpts',modelOpts);

%% VALIDATION
% Note that the steps below will over-write the "stimulus" and "data"
% variables. We may wish to make this more robust later.

% Load the validation dataset
fileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),'data',[sub '_city1B_stimulus_data_bMask-city1AB-100.mat']);
load(fileName,'stimulus','data')

% Re-create the model, now with the new stimulus and data
model = heading(data,stimulus,tr,modelOpts{:});

% Prepare the validation data into a signal. You would want to loop over
% vertices to test this across the entire dataset.
ii = 1; % Just do the first vertex
dataPrep = catcell(2,model.prep(data));
signal = model.clean(dataPrep(vxs(ii),:)');

% Obtain the parameter values of the previous fit for this vertex
x = results.params(vxs(ii),:);

% Obtain the R2 metric for the signal, given the params (x) we have
% previously calculated
valR2 = model.metric(signal,x);


