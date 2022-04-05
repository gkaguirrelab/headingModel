% simSignalDemo.m
%
% This script demonstrates creating a simulated signal for a model, adding
% some noise, and testing if we can recover the model parameters
%

% Load some example data and stimulus file to base the simulation upon
sub='sub-08';
fileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),'data',[sub '_city1A_stimulus_data_bMask-city1AB-100.mat']);

% Load the stimulus and data variables
load(fileName,'stimulus','data')


%% Simulated a heading direction
% This code replaces the actual stimulus with a fully randomized ordering
% of heading values, sampled from a specified number of discrete headings.
nTRs = size(stimulus{1},2); % TRs per acquisition
simBins = 8;
binSeparation = (2*pi/simBins);
binCenters = 0:binSeparation:(2*pi)-binSeparation;
for ii=1:length(stimulus)
    stimulus{ii} = datasample(binCenters,nTRs);
end


% The TR of the experiment, in seconds
tr=2;

% Define modelOpts
nFilterBins = 8;
modelOpts = {'nFilterBins',nFilterBins};

% Create a model object
model = heading(data,stimulus,tr,modelOpts{:});

% Get the x0 param values
x0 = model.initial;

% Set the gain of the adaptation effect to zero.
x0(1) = 0;

% Pick a preferred direction, which is the bin center closest to the
% preferred direction
preferredDirection = pi/2;
binSeparation = (2*pi/nFilterBins);
binCenters = 0:binSeparation:(2*pi)-binSeparation;
[~,idx]=min(abs(binCenters - preferredDirection));
preferredDirection = binCenters(idx);
fprintf('The preferred direction is: %2.2f rads\n',preferredDirection)

% Set the parameter for the preferred direction to unity
[~,idx] = min(abs(binCenters - preferredDirection));
myBin = 4+idx; % The first four parameters handle the adaptation model
x0(myBin) = 1;

% Get the simulated signal for the x params
simSignal = model.forward(x0);

% Plot this
figure
subplot(4,1,1);
thisVec = stimulus{1};
plot(thisVec,'-k');
hold on
idx = find(thisVec==preferredDirection);
plot(idx,preferredDirection,'*r')
xlabel('time [TRs]');
ylabel('heading [rads]');
title(sprintf('simulated heading direction (bins = %d)',simBins));
set(gca,'YTick',0:pi/2:2*pi) 
set(gca,'YTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'})

subplot(4,1,2);
plot(simSignal(1:nTRs));
xlabel('time [TRs]');
ylabel('BOLD response');
title(sprintf('simulated BOLD response (bins = %d)',nFilterBins));

% Let's see if we can recover these parameters. First, create a data
% variable from the simSignal
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
plot(binCenters,x0(5:5+nFilterBins-1),'-k');
hold on
plot(binCenters,x1(5:5+nFilterBins-1),'*r');
set(gca,'XTick',0:pi/2:2*pi) 
set(gca,'XTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'})
ylabel('model parameter');
xlabel('bin centers [rads]');
title('simulated and recovered bin amplitude');

