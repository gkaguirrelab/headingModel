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

%% Optional
% This code replaces the actual stimulus with a simplified version. The
% heading direction is set to be pi/2 a, in
% which the heading direction rotates through the available angles over a
% cycle time of 45 TRs.
simBins = 45;
preferredDirection = pi/2;
binSeparation = (2*pi/simBins);
binCenters = 0:binSeparation:(2*pi)-binSeparation;
[~,idx]=min(abs(binCenters - preferredDirection));
preferredDirection = binCenters(idx);
nonPreferredDirections = binCenters((1:simBins)~=idx);
thisVec = nan(1,330);
idx = logical(zeros(1,330));
idx(1:simBins:330) = 1;
thisVec(idx) = preferredDirection;
thisVec(~idx) = datasample(nonPreferredDirections,sum(~idx));
for ii=1:length(stimulus)
    stimulus{ii} = thisVec;
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

% We are working in % change units, so set the first (gain) parameter of
% the adaptation effect to a value in units of % change (something between
% 0 and 0.5 would be typical)
x0(1) = 0;

% The filter bin parameters are at index positions [5:5+nFilterBins-1]. Set
% a value for one of the bins
binSeparation = (2*pi/nFilterBins);
binCenters = 0:binSeparation:(2*pi)-binSeparation;
[~,idx] = min(abs(binCenters - preferredDirection));
myBin = 4+idx;
x0(myBin) = 1;

% Get the simulated signal for the x params
simSignal = model.forward(x0);

% Plot this
figure
subplot(3,1,1);
plot(thisVec,'-k');
hold on
idx = find(thisVec==preferredDirection);
plot(idx,preferredDirection,'*r')
xlabel('time [TRs]');
ylabel('heading [rads]');
title(sprintf('simulated heading direction (bins = %d)',simBins));
set(gca,'YTick',0:pi/2:2*pi) 
set(gca,'YTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'})

subplot(3,1,2);
plot(simSignal(1:330));
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
subplot(3,1,3);
plot(fitSignal(1:330));
xlabel('time [TRs]');
ylabel('BOLD response');
title(sprintf('fitted BOLD response (bins = %d)',nFilterBins));

% Compare the simulated and recovered values
fprintf('Overall gain [simulated vs recovered]: %2.2f vs. %2.2f \n',x0(1),x1(1));
fprintf('Overall exponent [simulated vs recovered]: %2.2f vs. %2.2f \n',x0(2),x1(2));
fprintf('Overall time constant [simulated vs recovered]: %2.2f vs. %2.2f \n',x0(4),x1(4));
fprintf('Heading bins: \n');
fprintf('  simulated: %2.1f, %2.1f, %2.1f, %2.1f, %2.1f, %2.1f, %2.1f, %2.1f \n',x0(5:5+nFilterBins-1));
fprintf('  recovered: %2.1f, %2.1f, %2.1f, %2.1f, %2.1f, %2.1f, %2.1f, %2.1f \n',x1(5:5+nFilterBins-1));
fprintf('HRF params: \n');
fprintf('  simulated: %2.2f, %2.2f, %2.2f \n',x0(5+nFilterBins:5+nFilterBins+2));
fprintf('  recovered: %2.2f, %2.2f, %2.2f \n',x1(5+nFilterBins:5+nFilterBins+2));

% Could try adding some noise to simSignal and see how well you recover the
% parameters
