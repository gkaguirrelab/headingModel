% simSignalUniWts.m
% tbUseProject('headingModel')
% addpath(genpath('/Users/zhenganglu/Documents/MATLAB/toolboxes'));
% This script demonstrates creating a simulated signal for a model, adding
% some noise, and testing if we can recover the model parameters
%
% Flags that control simulation properties
useRealHeadingFlag = true;
uniSimulationFlag = true;
hrfSearchFlag = true;

% Fixed variables of the simulation
tr = 2;             % The TR of the experiment, in seconds
preferredDirection = pi/2;

simSigma = pi/48;      % The width of the distribution of the bin weights when
                    %    simulating a distribution of heading weights
simBins = 45;       % how many bins to simulate in the signal generation
nFilterBins = 45;   % how many filters in the decoding model
% FWHM = 2*pi/simBins; % width of filter in term of degree/radian
% kappa = 8*log(2)/FWHM^2; % convert to width
nFixedParams = 2; % corresponding to the adaptation gain and epsilon
% Pick an example subject on which to base the demo
sub='sub-08';

% Load the stimulus and data variables
fileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),['data/' sub] ...
    ,[sub '_city1A_stimulus_data_bMask-100.mat']);
load(fileName,'stimulus','data')

%% Simulate a heading direction
% This code replaces the actual stimulus with a fully randomized ordering
% of heading values, sampled from a specified number of discrete headings,
% under the control of a flag
nTRs = size(stimulus{1},2); % TRs per acquisition
simBinSeparation = (2*pi/simBins);
simBinCenters = 0:simBinSeparation:(2*pi)-simBinSeparation;

if ~useRealHeadingFlag
    for ii=1:length(stimulus)
        stimulus{ii} = datasample(simBinCenters,nTRs);
    end
end

% Define modelOpts for the simulation model

modelOpts = {'nFilterBins',simBins,'hrfSearch',false,'typicalGain',1};

% Create a model object
modelSim = heading(data,stimulus,tr,modelOpts{:});

% Get the x0 param values
x0 = modelSim.initial;

% Set the gain and epsilon of the adaptation effect to some values
x0(1:2) = [1 0.8];

% Two choices for modeling the heading direction effect: "one hot", or a
% circular Gaussian distribution of amplitudes, under the control of a flag
if uniSimulationFlag

    % all the heading bins have the same weights
    myBin = nFixedParams+nFilterBins; % The first four parameters handle the adaptation model
    x0(nFixedParams+1:myBin) = 1;
else
    % Create a von Misses distribution of bin weights in the circular
    % heading space. The function circ_vmpdf takes the bin centers, a
    % preferred direction, and the "concentration" parameter kappa, which
    % has an inverse relationship to the sigma of a conventional Gaussian
    % distribution
    kappa = sqrt(1/simSigma^2);
    x0(nFixedParams+1:simBins+nFixedParams) = circ_vmpdf(simBinCenters,preferredDirection,kappa);
end

% Get the simulated signal for the x params
simSignal = modelSim.forward(x0);
%% Add noise to the signal
simSignalNoise = simSignal+0.5*randn(size(simSignal));

%% Let's see if we can recover these parameters
% First, create a data variable from the simSignal
nAcq = size(data,2);
nTRsPerAcq = length(simSignal)/nAcq;
for ii=1:nAcq
    data{ii} = simSignalNoise(nTRsPerAcq*(ii-1)+1:nTRsPerAcq*ii)';
end

% Update the modelOpts with the number of bins used for decoding
modelOpts = {'nFilterBins',nFilterBins,'hrfSearch',hrfSearchFlag,'typicalGain',1};

% Call the forwardModel
results = forwardModel(data,stimulus,tr,'modelClass','heading','vxs',1,'modelOpts',modelOpts);

% Pull out the params
x1 = results.params;

% Obtain the model fit
modelOut = heading(data,stimulus,tr,modelOpts{:});
fitSignal = modelOut.forward(x1);

%%
% Plot this
xWts=x1(1,nFixedParams+1:simBins+nFixedParams);
figure
subplot(4,1,1);
thisVec = stimulus{1};
plot(thisVec,'-k');
hold on
xlabel('time [TRs]');
ylabel('heading [rads]');
title('real heading direction');
set(gca,'YTick',0:pi/2:2*pi) 
set(gca,'YTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'})

subplot(4,1,2);
plot(simSignalNoise(1:nTRs));
xlabel('time [TRs]');
ylabel('BOLD response');
title(sprintf('simulated BOLD response (bins = %d)',simBins));

% Add this to the plot
subplot(4,1,3);
plot(fitSignal(1:nTRs));
xlabel('time [TRs]');
ylabel('BOLD response');
title(sprintf('fitted BOLD response (bins = %d)',nFilterBins));

% Compare the simulated and recovered values. First create the bin centers
% used for model read-out
modelBinSeparation = (2*pi/nFilterBins);
modelBinCenters = 0:modelBinSeparation:(2*pi)-modelBinSeparation;

subplot(4,1,4);
plot(simBinCenters,x0(nFixedParams+1:nFixedParams+simBins),'-k');
hold on
plot(modelBinCenters,xWts,'*r');
% plot(modelBinCenters,x1(nFixedParams+1:nFixedParams+nFilterBins),'*r');
set(gca,'XTick',0:pi/2:2*pi) 
set(gca,'XTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'})
ylim([1-0.3, 1+0.3]);
ylabel('model parameter');
xlabel('bin centers [rads]');
title('simulated and recovered bin amplitude');
