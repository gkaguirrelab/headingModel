% load stimulus and fMRI data
% tbUseProject('headingModel')

fileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),'data','sub-00_stimulus_data.mat');
load(fileName,'stimulus','data')
tr=2;

%% Call the forwardModel
results = forwardModel(data,stimulus,tr,'modelClass','heading');

% Show the results figures
figFields = fieldnames(results.figures);
if ~isempty(figFields)
    for ii = 1:length(figFields)
        figHandle = struct2handle(results.figures.(figFields{ii}).hgS_070000,0,'convert');
        set(figHandle,'visible','on')
    end
end