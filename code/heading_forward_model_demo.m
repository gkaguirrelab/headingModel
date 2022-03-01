% load stimulus and fMRI data
% tbUseProject('headingModel')

% The mat file we need is located in the "data" directory. Construct a path
% to this file, starting with the location of this script.
sub='sub-08';
fileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),'data',[sub '_stimulus_data.mat']);

% Load the stimulus and data variables
load(fileName,'stimulus','data')

% The TR of the experiment, in seconds
tr=2;

% Define a vxs if we wish to analyze just a subset of available vertices.
% By default, analyze them all
vxs = 1:size(data{1},1);

% Call the forwardModel
results = forwardModel(data,stimulus,tr,'modelClass','heading','vxs',vxs);

% Show the results figures
figFields = fieldnames(results.figures);
if ~isempty(figFields)
    for ii = 1:length(figFields)
        figHandle = struct2handle(results.figures.(figFields{ii}).hgS_070000,0,'convert');
        set(figHandle,'visible','on')
    end
end