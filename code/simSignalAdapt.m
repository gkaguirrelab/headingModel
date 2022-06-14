% simSignalAdapt.m
% tbUseProject('headingModel')
% This script demonstrates creating a simulated signal for a model, adding
% some noise, and testing if we can recover the model parameters
%
clear;
close all;
% Flags that control simulation properties
useRealHeadingFlag = true;
oneHotSimulationFlag = false;

% Load some example data and stimulus file to base the simulation upon
sub='sub-08';

% Load the stimulus and data variables
fileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),['data/' sub] ...
    ,[sub '_city1A_stimulus_data_bMask-100.mat']);
load(fileName,'stimulus','data')
% 
%% Simulated a heading direction
% This code replaces the actual stimulus with a fully randomized ordering
% of heading values, sampled from a specified number of discrete headings,
% under the control of a flag
nTRs = size(stimulus{1},2); % TRs per acquisition
if ~useRealHeadingFlag
    simBins = 8; % how many unique direction are there
    binSeparation = (2*pi/simBins);
    binCenters = 0:binSeparation:(2*pi)-binSeparation;
    for ii=1:length(stimulus)
        stimulus{ii} = datasample(binCenters,nTRs);
    end
end

%% Create a model object
tr = 2; % The TR of the experiment, in seconds
% Define modelOpts
nFixedParams = 2; % corresponding to the adaptation gain and epsilon
nFilterBins = 16; % how many filters in the model
modelOpts = {'nFilterBins',nFilterBins,'hrfSearch',false,'typicalGain',1};

binSeparation = (2*pi/nFilterBins);
binCenters = 0:binSeparation:(2*pi)-binSeparation;
% Pick a preferred bin, which is the bin center closest to the preferred
% direction
% store which bin in the simulated heading direction is closest to the preferred direction
% Set a preferredDirection
preferredDirection = pi/2; %2*binSeparation;%pi/2; %7*pi/4; %pi/2;

model = heading(data,stimulus,tr,modelOpts{:});
% Get the x0 param values
x0 = model.initial;
% Set the gain of the adaptation effect to zero, or 25% of heading change
x0(1) = 0;
% x0(1:2) = [1 0.5];
% 
% Two choices for modeling the heading direction effect: "one hot", or a
% circular Gaussian distribution of amplitudes, under the control of a flag
if oneHotSimulationFlag
    [~,idx]=min(abs(binCenters - preferredDirection));
    preferredDirection = binCenters(idx);
    preferredDirectionInHeadingVector = binCenters(idx);
    myBin = nFixedParams+idx; % The first four parameters handle the adaptation model
    x0(myBin) = 1;
else
    x0(nFixedParams+1:nFilterBins+nFixedParams) = circ_vmpdf(binCenters,preferredDirection,binSeparation);
end
fprintf('The preferred direction is: %2.2f degree\n', ...
    180*preferredDirection/pi)
% Get the simulated signal for the x params
simSignal = model.forward(x0);

% Plot this
figure
subplot(4,1,1);
thisVec = stimulus{1};
plot(thisVec,'-k');
hold on
idx = find(round(thisVec,2)==round(preferredDirection,2));
% idx = find(round(thisVec,2)==round(preferredDirectionInHeadingVector,2));
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

%% Let's see if we can recover these parameters. First, create a data
% variable from the simSignal
nAcq = size(data,2);
nTRsPerAcq = length(simSignal)/nAcq;
simBOLD=[];
for ii=1:nAcq
    data{ii} = simSignal(nTRsPerAcq*(ii-1)+1:nTRsPerAcq*ii)';
    simBOLD(:,ii)=simSignal(nTRsPerAcq*(ii-1)+1:nTRsPerAcq*ii)';
end
%%
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
%% Try to recover with different numbers of heading bins
FilterWidth = [8 10 15 20 22.5 24 30 36 45 60];
nFilterBinsTests = 360./FilterWidth;
X={};
for i=1:length(nFilterBinsTests)
    modelOptsTest = {'nFilterBins',nFilterBinsTests(i),'hrfSearch',false,'typicalGain',1};
    % Call the forwardModel
    results = forwardModel(data,stimulus,tr,'modelClass','heading','vxs',1,'modelOpts',modelOptsTest);

    % Pull out the params
%     X{i}=results.params(nFixedParams+1:nFixedParams+nFilterBins);
    X{i}=results.params;
end

%% save simulated data for linear encoding model
simStimuli = stimulus;
outPath='/Users/zhenganglu/Projects/brainSLAM/fMRI_DATA/sub-sim';
fn=[outPath '/sub-08_sim-adapt-' num2str(x0(1)) '-' num2str(x0(2))...
    '_head-' num2str(nFilterBins)];
figName=[fn '_gaussian_stim-real.png']
dataName = [fn '_gaussian_stim-real.mat']
if oneHotSimulationFlag
    figName=[fn '_one-hot_stim-real.png'];
    dataName = [fn '_one-hot_stim-real.mat'];
end
saveas(gcf,figName);
save(dataName, 'simSignal', 'simBOLD', 'simStimuli', 'x0', 'X', 'fitSignal');
