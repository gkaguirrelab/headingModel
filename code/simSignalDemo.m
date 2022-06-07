% simSignalDemo.m
% tbUseProject('headingModel')
% This script demonstrates creating a simulated signal for a model, adding
% some noise, and testing if we can recover the model parameters
%

% Flags that control simulation properties
useRealHeadingFlag = false;
oneHotSimulationFlag = true;

% Pick an example subject on which to base the demo
sub='sub-08';

% Load the stimulus and data variables
fileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),['data/' sub] ...
    ,[sub '_city1A_stimulus_data_bMask-city1AB-100.mat']);
load(fileName,'stimulus','data')


%% Simulate a heading direction

% This code replaces the actual stimulus with a fully randomized ordering
% of heading values, sampled from a specified number of discrete headings,
% under the control of a flag
nTRs = size(stimulus{1},2); % TRs per acquisition
simBins = 16; % how many unique direction are there
binSeparation = (2*pi/simBins);
binCenters = 0:binSeparation:(2*pi)-binSeparation;
if ~useRealHeadingFlag
    for ii=1:length(stimulus)
        stimulus{ii} = datasample(binCenters,nTRs);
    end
end

% Set a preferredDirection, and store which bin in the simulated heading
% direction is closest to the preferred direction
preferredDirection = pi/2; %2*binSeparation;%pi/2; %7*pi/4; %pi/2;
[~,idx]=min(abs(binCenters - preferredDirection));
preferredDirectionInHeadingVector = binCenters(idx);

% The TR of the experiment, in seconds
tr = 2;

% Define modelOpts
nFixedParams = 2; % corresponding to the adaptation gain and epsilon
nFilterBins = 16; % how many filters in the model
filterWidth=360/nFilterBins;
modelOpts = {'nFilterBins',nFilterBins,'hrfSearch',false,'typicalGain',1};

% Create a model object
model = heading(data,stimulus,tr,modelOpts{:});

% Get the x0 param values
x0 = model.initial;

% Set the gain and epsiln of the adaptation effect to some values
x0(1:2) = [0.5 1.2];

% Two choices for modeling the heading direction effect: "one hot", or a
% circular Gaussian distribution of amplitudes, under the control of a flag
if oneHotSimulationFlag
    % Pick a preferred bin, which is the bin center closest to the preferred
    % direction
    binSeparation = (2*pi/nFilterBins);
    binCenters = 0:binSeparation:(2*pi)-binSeparation;
    [~,idx]=min(abs(binCenters - preferredDirection));
    preferredDirection = binCenters(idx);
    fprintf('The preferred direction is: %2.2f degree\n', ...
        180*preferredDirection/pi)
    myBin = nFixedParams+idx; % The first four parameters handle the adaptation model
    x0(myBin) = 1;
else
    x0(nFixedParams+1:nFilterBins+nFixedParams) = circ_vmpdf(binCenters,preferredDirection,binSeparation);
end

% Get the simulated signal for the x params
simSignal = model.forward(x0);

% Plot this
figure
subplot(4,1,1);
thisVec = stimulus{1};
plot(thisVec,'-k');
hold on
idx = find(round(thisVec,2)==round(preferredDirectionInHeadingVector,2));
plot(idx,preferredDirection,'*r')
xlabel('time [TRs]');
ylabel('heading [rads]');
% title(sprintf('simulated heading direction (bins = %d)',simBins));
title('real heading direction');
set(gca,'YTick',0:pi/2:2*pi) 
set(gca,'YTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'})

subplot(4,1,2);
plot(simSignal(1:nTRs));
xlabel('time [TRs]');
ylabel('BOLD response');
title(sprintf('simulated BOLD response (bins = %d)',nFilterBins));

%% Let's see if we can recover these parameters
% First, create a data variable from the simSignal
nAcq = size(data,2);
nTRsPerAcq = length(simSignal)/nAcq;
for ii=1:nAcq
    data{ii} = simSignal(nTRsPerAcq*(ii-1)+1:nTRsPerAcq*ii)';
end

% Call the forwardModel
results = forwardModel(data,stimulus,tr,'modelClass','heading','vxs',1,'modelOpts',modelOpts);

% Pull out the params
x1 = results.params;

% Obtain the model fit
fitSignal = model.forward(x1);

% Add this to the plot
subplot(4,1,3);
plot(fitSignal(1:nTRs));
xlabel('time [TRs]');
ylabel('BOLD response');
title(sprintf('fitted BOLD response (bins = %d)',nFilterBins));

% Compare the simulated and recovered values
subplot(4,1,4);
plot(binCenters,x0(nFixedParams+1:nFixedParams+nFilterBins),'-k');
hold on
plot(binCenters,x1(nFixedParams+1:nFixedParams+nFilterBins),'*r');
set(gca,'XTick',0:pi/2:2*pi) 
set(gca,'XTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'})
ylabel('model parameter');
xlabel('bin centers [rads]');
title('simulated and recovered bin amplitude');
