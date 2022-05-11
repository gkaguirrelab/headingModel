% simSignalDemo.m
% tbUseProject('headingModel')
% This script demonstrates creating a simulated signal for a model, adding
% some noise, and testing if we can recover the model parameters
%

% Load some example data and stimulus file to base the simulation upon
sub='sub-08';
% The TR of the experiment, in seconds
tr=2;
% Define modelOpts
nFilterBins = 36; % how many filters in the model
filterWidth=360/nFilterBins;
modelOpts = {'nFilterBins',nFilterBins};
one_hot=1;
realHD=1;
% Load the stimulus and data variables
fileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),['data/' sub] ...
    ,[sub '_city1A_stimulus_data_bMask-city1AB-100.mat']);
load(fileName,'stimulus','data')
% 
%% Simulated a heading direction
% This code replaces the actual stimulus with a fully randomized ordering
% of heading values, sampled from a specified number of discrete headings.
nTRs = size(stimulus{1},2); % TRs per acquisition
simBins = 8; % how many unique direction are there
binSeparation = (2*pi/simBins);
binCenters = 0:binSeparation:(2*pi)-binSeparation;
if realHD==0
    for ii=1:length(stimulus)
        stimulus{ii} = datasample(binCenters,nTRs);
    end
end

% Set a preferredDirection, and store which bin in the simulated heading
% direction is closest to the preferred direction
preferredDirection = pi/2; %2*binSeparation;%pi/2; %7*pi/4; %pi/2;
[~,idx]=min(abs(binCenters - preferredDirection));
preferredDirectionInHeadingVector = binCenters(idx);

% Create a model object
model = heading(data,stimulus,tr,modelOpts{:});

% Get the x0 param values
x0 = model.initial;

% Set the gain of the adaptation effect to zero.
x0(1) = 0;
% 
% % Pick a preferred direction, which is the bin center closest to the
% % preferred direction
binSeparation = (2*pi/nFilterBins);
binCenters = 0:binSeparation:(2*pi)-binSeparation;
[~,idx]=min(abs(binCenters - preferredDirection));
preferredDirection = binCenters(idx);
fprintf('The preferred direction is: %2.2f degree\n', ...
    180*preferredDirection/pi)

% Set the parameter for the preferred direction to unity, all the other to
% have values following circular gaussion profile
% one hot encoding:
% one_hot=[0,0,1,0,0,0,0,0];imagesc(one_hot);colormap(jet);axis off;colorbar

[~,idx] = min(abs(binCenters - preferredDirection));
myBin = 4+idx; % The first four parameters handle the adaptation model
x0(myBin) = 1;
% gaussian encoding: gau=circ_vmpdf(binCenters,preferredDirection,binSeparation); imagesc(gau');colormap(jet)
if one_hot==0
    x0(5:nFilterBins-1+5) = circ_vmpdf(binCenters,preferredDirection,binSeparation);
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

%% Let's see if we can recover these parameters. First, create a data
% variable from the simSignal
nAcq = size(data,2);
nTRsPerAcq = length(simSignal)/nAcq;
simBOLD=[];
for ii=1:nAcq
    data{ii} = simSignal(nTRsPerAcq*(ii-1)+1:nTRsPerAcq*ii)';
    simBOLD(:,ii)=simSignal(nTRsPerAcq*(ii-1)+1:nTRsPerAcq*ii)';
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
%% Try to recover with different numbers of heading bins
FilterWidth = [8 10 15 20 24 30 36 45 60];
nFilterBinsTests = 360./FilterWidth;
X={};
for i=1:length(nFilterBinsTests)
    
    modelOptsTest = {'nFilterBins', nFilterBinsTests(i)};
    % Call the forwardModel
    results = forwardModel(data,stimulus,tr,'modelClass','heading','vxs',1,'modelOpts',modelOptsTest);

    % Pull out the params
    X{i}=results.params(5:nFilterBinsTests(i)-1+5);
end

% Obtain the model fit
fitSignal = model.forward(x1);
% save([outPath '/sub-08_head-8-45_simSignalStim-' num2str(nSim) '_circ.mat'], ...
%     'simSignal', 'simBOLD', 'simStimuli', 'x0', 'x1', 'fitSignal');
%% save simulated data for linear encoding model
simStimuli = stimulus;
outPath='/Users/zhenganglu/Projects/brainSLAM/fMRI_DATA/sub-sim';
fn=[outPath '/sub-08_sim_head-' num2str(nFilterBins) '-' num2str(filterWidth)];
figName=[fn '_gaussian_stim-real.png'];
dataName = [fn '_gaussian_stim-real.mat'];
if one_hot==1
    figName=[fn '_one-hot_stim-real.png'];
    dataName = [fn '_one-hot_stim-real.mat'];
end
saveas(gcf,figName);
save(dataName, 'simSignal', 'simBOLD', 'simStimuli', 'x0', 'X', 'fitSignal');
