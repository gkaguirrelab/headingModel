% load stimulus and fMRI data
% tbUseProject('headingModel')

% The mat file we need is located in the "data" directory. Construct a path
% to this file, starting with the location of this script.
fileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),'data','sub-00_stimulus_data.mat');

% Load the stimulus and data variables
load(fileName,'stimulus','data')

% The TR of the experiment, in seconds
tr=2;

% Call the forwardModel
results = forwardModel(data,stimulus,tr,'modelClass','heading');

% Show the results figures
figFields = fieldnames(results.figures);
if ~isempty(figFields)
    for ii = 1:length(figFields)
        figHandle = struct2handle(results.figures.(figFields{ii}).hgS_070000,0,'convert');
        set(figHandle,'visible','on')
    end
end