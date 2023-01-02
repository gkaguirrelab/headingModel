function [x0, x1] = simEngine(noiseScale,binWeightMax,fixedParamVector,lassoRegularization,options)
% Function to conduct simulations using the heading forwardModel
%
% Syntax:
%   [x0, x1] = simEngine(noiseScale,binWeightMax,adaptationGain,lassoRegularization,options)
%
% Description:
%   Supports various simulations using the heading forwardModel. See
%   examples below.
%
% Inputs:
%   noiseScale            - The amount of white noise to add to the
%                           simulated data
%   binWeightMax          - The max bin weight for the simulated response
%                           to heading direction. Set to zero to have
%                           no response to heading direction.
%   fixedParamVector      - A vector of the values to use for the fixed
%                           parameters in the simulation, with the entries
%                           corresponding to:
%                             - adaptation gain
%                             - adaptation exponent
%                             - adaptation tau
%   lassoRegularization   - The regularization on the lasso penalty that is
%                           derived from the filter bin weights
%
% Optional key/value pairs:
%   nSimBins              - The number of filter bins to use in creating
%                           the simulated data
%   nFitBins              - The number of filter bins to use in fitting
%                           the simulated data
%   padTimeSecs           - The number of seconds to pad the start of the 
%                           stimulus sequence, to provide some "warm up"
%                           for the model
%   headingDownsampleFactor - Scalar index. How much down-sampling to
%                           perform on the heading direction vector.
%  'useRealHeading'       - Flag to control using the actual subject
%                           heading vectors. If set to false a permuted
%                           version of the actual heading is used.
%  'hrfSearch'            - Flag to control if the fitting is performed
%                           with a search over the parameters of the hrf
%  'adaptSearch'          - Flag to control if the model implements an
%                           adaptation effect.
%  'makePlots'            - Flag to control if plots are produced
%
% Outputs:
%   x0, x1            - Vectors. The modeled and recovered params.
%
% Examples:
%{
    % No heading response, no adaptation, with some small amount
    % of noise and no lasso penalty. This causes the bin weights to show
    % an oscillation
    fixedParams = [0 1];
    simEngine(1,0,fixedParams,0);
%}
%{
    % No heading response, no adaptation, with large noise, but include the
    % lasso penalty. This removes the oscillation in the bin weights. Note
    % that y-axis scale for the plot of bin weights is quite small.
    fixedParams = [0 1];
    simEngine(4,0,fixedParams,0.05);
%}
%{
    % Now demonstrate that with a veridical heading effect, the fit
    % recovers the bin weights well, even in the presence of substantial
    % noise
    fixedParams = [0 1 5];
    simEngine(4,1,fixedParams,0.05);
%}
%{
    % Recover adaptation and heading direction effects
    fixedParams = [1 1 15];
    [x0,x1]=simEngine(4,0.5,fixedParams,0.05);
%}
%{
    % Recover adaptation and heading direction effects for the downsampled 
    % heading direction vector
    fixedParams = [1 1 15];
    [x0,x1]=simEngine(4,0.5,fixedParams,0.05,'headingDownsampleFactor',3);
%}
%{
    % Set some fixed parameters, and see if we can recover them, including
    % an interpolated preferred heading direction.
    fixedParams = [1.5 1 15];
    nBins = 45;
    [x0, x1] = simEngine(1,1,fixedParams,'hrfSearch',false);
    fprintf('simulated and recovered adaptation gain: [%2.2f, %2.2f] \n',x0(1),x1(1));
    fprintf('simulated and recovered adaptation exponent: [%2.2f, %2.2f] \n',x0(2),x1(2));
    fprintf('simulated and recovered tau: [%2.2f, %2.2f] \n',x0(3),x1(3));
    % Obtain the interpolated peak of the preferred heading direction
    binSupportFit = 1:0.01:nBins;
    directionSupportFit = linspace(0,2*pi-(2*pi/nBins),length(binSupportFit));
    x0Fit = spline(1:nBins,x0(4:end-3),binSupportFit);
    x1Fit = spline(1:nBins,x1(4:end-3),binSupportFit);
    % Report the peak preferred bin
    [~,idx0]=max(x0Fit);
    fprintf('The interpolated, veridical heading direction is %2.2f \n',directionSupportFit(idx0));
    [~,idx1]=max(x1Fit);
    fprintf('The interpolated, recovered heading direction is %2.2f \n',directionSupportFit(idx1));
    signedError = directionSupportFit(idx0)-directionSupportFit(idx1);
    fprintf('The signed error in recovered heading direction is %2.2f \n',signedError);    
%}


arguments
    noiseScale {isscalar,mustBeNumeric} = 1
    binWeightMax {isscalar,mustBeNumeric} = 1
    fixedParamVector {isvector,mustBeNumeric} = [0, 1, 0]
    lassoRegularization {isscalar,mustBeNumeric} = 0.05
    options.nSimBins {isscalar,mustBeNumeric} = 45;       % how many bins to simulate in the signal generation
    options.nFitBins {isscalar,mustBeNumeric} = 45;   % how many filters in the decoding model
    options.padTimeSecs {isscalar,mustBeNumeric} = 16;   % how many filters in the decoding model
    options.headingDownsampleFactor = 1;
    options.useRealHeading {islogical} = true
    options.hrfSearch {islogical} = true
    options.adaptSearch {islogical} = true
    options.makePlots {islogical} = true
end


%% Fixed variables of the simulation
%   - tr in seconds
%   - The preferred direction to be modeled, in radians
%   - The number of TRs to extend the first heading value prior to the
%     start ofthe scan. This is necessary to avoid a ramp up of BOLD fMRI
%     response at the start of the experiment. This could also be
%     accomplished by discarding the first n TRs of the data in the
%     modeling effort.
% Data TRs per acquisition
tr = 2;
nTRs = 330;
preferredDirection = pi;
nFixedParams = 3; % corresponding to the adaptation gain, epsilon, and muu

%% Create a stimulus
% We load a set of heading direction vectors from an example subject, and
% either use these values, or base a simulation upon them.
fileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),'data/', ...
    'sub-08_city1A_stimulus.mat');
load(fileName,'stimulus')
nDTs = size(stimulus{1},2);
deltaTStim = 660/3384; % The temporal sampling interval (in seconds) of the
                       % heading direction vector

% Resample the heading direction vector if requested
if options.headingDownsampleFactor ~= 1
    stimulus = cellfun(@(x) wrapTo2Pi(decimate(unwrap(x),options.headingDownsampleFactor)),stimulus,'UniformOutput',false);
    nDTs = length(stimulus{1});
    deltaTStim = deltaTStim * options.headingDownsampleFactor;
end

% How many samples do we need to pad the stimulus to have the requested
% duration of padding in seconds?
nDTsToPad = round(options.padTimeSecs*(1/deltaTStim));

% Either use the actual heading direction, or randomize the order of the
% heading values in each acquisition
if ~options.useRealHeading
    for ii=1:length(stimulus)
        stimulus{ii} = stimulus{ii}(randperm(nTRs));
    end
end

% Extend the stimulus back in time to provide some "warm up" for the HRF;
% create a stimTime array
stimTime={};
for ii=1:length(stimulus)

    % Grab an acquisition from the cell array
    thisStim = stimulus{ii};

    % Create a stimTime vector
    stimTime{ii} = -nDTsToPad*deltaTStim:deltaTStim:deltaTStim*(length(thisStim)-1);

    % Pad this stim
    thisStim = [repmat(thisStim(1),1,nDTsToPad), thisStim];
    stimulus{ii} = thisStim;

end

%% Create data

% Define the dataTime
dataTime = tr*(0:nTRs-1);
headingTime= -nDTsToPad*deltaTStim:deltaTStim:deltaTStim*(nDTs-1);

% Define modelOpts for the simulation model
modelOpts = {'nFilterBins',options.nSimBins,'typicalGain',1};

% To initialize the model, we need to pass some "dummyData" that informs
% the model of the dimensions of the fit to be returned
dummyData = repmat({zeros(1, nTRs)},1,size(stimulus,2));

% Create a model object
modelSim = heading(dummyData,stimulus,tr,'stimTime',stimTime,modelOpts{:});

% Get the x0 param values
x0 = modelSim.initial;

% What's the total number of params?
nParams = length(x0);

% Set the fixed params to the passed values
x0(1:length(fixedParamVector)) = fixedParamVector;

% Create a von Misses distribution of bin weights in the circular heading
% space. The function circ_vmpdf takes the bin centers, a preferred
% direction, and the "concentration" parameter kappa, which has an inverse
% relationship to the sigma of a conventional Gaussian
% distribution
simBinSeparation = (2*pi/options.nSimBins);
simBinCenters = 0:simBinSeparation:(2*pi)-simBinSeparation;
FWHM = simBinSeparation;
sigma = FWHM/(2*sqrt(2*log(2)));
kappa = 1/(10*sigma^2);
x0(nFixedParams+1:options.nSimBins+nFixedParams) = binWeightMax.*circ_vmpdf(simBinCenters,preferredDirection,kappa);

% Get the simulated neural signal for the x params. To do so, we pass zeros
% for the hrf params
xDelta = x0;
xDelta(nParams-2:nParams) = 0;
simSignalNeural = modelSim.forward(xDelta);

% Get the simulated BOLD signal for the x params
simSignal = modelSim.forward(x0);

% Add noise to the BOLD signal
simSignalNoise = simSignal+noiseScale*randn(size(simSignal));


%% Recover parameters

% Create a data variable from the simSignal
nAcq = size(stimulus,2);
nTRsPerAcq = length(simSignal)/nAcq;
for ii=1:nAcq
    data{ii} = simSignalNoise(nTRsPerAcq*(ii-1)+1:nTRsPerAcq*ii)';
end

% Update the modelOpts with the number of bins used for decoding
modelOpts = {'nFilterBins',options.nFitBins,'hrfSearch',options.hrfSearch,'adaptSearch',options.adaptSearch,'typicalGain',1,'lassoRegularization',lassoRegularization};

% Call the forwardModel
results = forwardModel(data,stimulus,tr,'stimTime',stimTime,...
    'modelClass','heading','vxs',1,...
    'modelOpts',modelOpts);

% Pull out the params
x1 = results.params;

% Obtain the model fit
modelOut = heading(data,stimulus,tr,'stimTime',stimTime,modelOpts{:});
fitSignal = modelOut.forward(x1);
xAdaptWts=x1(1,1:nFixedParams);

%% Make plots
if options.makePlots
    xWts=x1(1,nFixedParams+1:options.nSimBins+nFixedParams);
    figure
    subplot(4,1,1);
    thisVec = stimulus{1};
    plot(headingTime(nDTsToPad:end),thisVec(nDTsToPad:end),'-ok');
    hold on
    xlabel('time [seconds]');
    ylabel('heading [rads]');
    title('real heading direction');
    set(gca,'YTick',0:pi/2:2*pi)
    set(gca,'YTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'})
    xlim([min(headingTime) max(headingTime)]);

    subplot(4,1,2);
    plot(dataTime,simSignalNeural(1:nTRs),'-ok');
    xlabel('time [seconds]');
    ylabel('Neural response');
    title(sprintf('simulated neural response downsampled to TRs (bins=%d, gain=%s, epsilon=%s, tau=%s)'...
        ,options.nSimBins,num2str(x0(1)), num2str(x0(2)), num2str(x0(3))));
    xlim([min(dataTime) max(dataTime)]);
    ylim([-1 2.5]);
    
    subplot(4,1,3);
    plot(dataTime,simSignalNoise(1:nTRs),'ok');
    hold on
    plot(dataTime,fitSignal(1:nTRs),'-or');
    xlabel('time [seconds]');
    ylabel('BOLD response');
    title(sprintf('simulated and fitted BOLD response (bins = %d)',options.nSimBins));
    xlim([min(dataTime) max(dataTime)]);
    ylim([-3 3]);

    % Create the bin centers used for model read-out
    modelBinSeparation = (2*pi/options.nFitBins);
    modelBinCenters = 0:modelBinSeparation:(2*pi)-modelBinSeparation;

    subplot(4,1,4);
    plot(simBinCenters,x0(nFixedParams+1:nFixedParams+options.nSimBins),'-ok');
    hold on
    plot(modelBinCenters,xWts,'or');
    set(gca,'XTick',0:pi/2:2*pi)
    set(gca,'XTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'})
    ylabel('model parameter');
    xlabel('bin centers [rads]');
    title('simulated and recovered bin amplitude');
    xlim([0 2*pi]);
end

end
