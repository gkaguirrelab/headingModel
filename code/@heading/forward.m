function [fit, hrf] = forward(obj, x)
% Forward model
%
% Syntax:
%   [fit, hrf] = obj.forward(x)
%
% Description:
%   Returns a time-series vector that is the predicted response to the
%   stimulus, based upon the parameters provided in x.
%
% Inputs:
%   x                     - 1xnParams vector.
%
% Optional key/value pairs:
%   none
%
% Outputs:
%   fit                   - 1xtime vector.
%   hrf                   - 1xn vector.
%
% zgl edit the kappa using findKappa function from vHD_encoding_model

% Obj variables
stimulus = obj.stimulus;
stimAcqGroups = obj.stimAcqGroups;
stimTime = obj.stimTime;
nParams = obj.nParams;
dataTime = obj.dataTime;
nFixedParamsAdapt = obj.nFixedParamsAdapt;
% nFixedParamsOther = obj.nFixedParamsOther;
filterResponse=obj.filterResponse;

% Break the parameters of x into named variables for code transparency
adaptGain = x(1);   % Gain of the adaptation effect
epsilon = x(2);     % Non-linear exponent of the neural signal
%tau = x(3);         % Time constant of the temporal integration
muu = x(3);         % proportion of prior updated to the current heading
%% Build a model of absolute heading direction
% We have loaded the stored filter bank from the object properties.
% Apply the filter weights as a single matrix multiplication.
FilterWeights=x(nFixedParamsAdapt+1:end-3);
neuralSignal = sum(filterResponse.*FilterWeights, 2);

% Obtain the temporal support for one the stimulus for one acquisition
stimTimeSingleAcq = stimTime(stimAcqGroups==1);

% Create an exponential kernel under the contol of tau (in units of
% seconds)
% exponentialIRF = exp(-1/tau*stimTimeSingleAcq);

% Normalize the kernel to have unit area
% exponentialIRF = exponentialIRF/sum(exponentialIRF);

% Define an empty heading change vector
headingChange = zeros(size(stimulus));

% Loop over acquisitions and obtain the absolute, circular value of the
% difference between the prior and heading aka adaptation effect. 
for run = 1:max(stimAcqGroups)
    temp = stimulus(stimAcqGroups==run);
    % headingChange(stimAcqGroups==run,:) = [0;abs(angdiff(temp))];
    % calculate the heading prior
    r0=temp(1);
    headingPrior = zeros(size(temp));
    headingPrior(1)=r0;
    for n=2:length(temp)
        currentHeading=temp(n);
        angChange=angdiff(r0,currentHeading);
        r=r0+(1-muu)*angChange;
        r0=r;
        headingPrior(n)=r;
    end
    % convert to 0-2pi
    headingPrior=wrapTo2Pi(headingPrior);
    % adaptation defined as the absolute ang difference between the prior 
    % and heading
    headingChange(stimAcqGroups==run) = abs(angdiff(headingPrior,temp));
end

% Apply an exponential parameter to produce a compressive or expansive
% non-linear mapping between heading change and neural response (a
% parameter value of 1 provides a linear mapping).
headingChange = headingChange.^epsilon;

% Apply the exponential kernel to the heading time series
% headingChange = conv2run(headingChange,exponentialIRF,stimAcqGroups);

% Scale the stimulus matrix by the gain parameter
neuralSignal = neuralSignal + headingChange*adaptGain;

% Create the HRF. First, check to see if we have all zeros for the hrf
% params, in which case the hrf is a delta function.
if all(x(nParams-2:nParams)==0)
    hrf = zeros(size(obj.flobsbasis,1),1);
    hrf(1)=1;
else
    switch obj.hrfType
        case 'flobs'
            hrf = makeFlobsHRF(x(nParams-2:nParams), obj.flobsbasis);
        case 'gamma'
            hrf = makeGammaHRF(x(nParams-2:nParams), obj.stimDeltaT);
        otherwise
            error('Not a recognized HRF type')
    end
end

% Store the hrf in the object
obj.hrf = hrf;

% Convolve the neuralSignal by the hrf, respecting acquisition boundaries
fit = conv2run(neuralSignal,hrf,stimAcqGroups);

% If the stimTime variable is not empty, resample the fit to match
% the temporal support of the data.
if ~isempty(stimTime)
    dataAcqGroups = obj.dataAcqGroups;
    dataTime = obj.dataTime;
    fit = resamp2run(fit,stimAcqGroups,stimTime,dataAcqGroups,dataTime);
end

% Apply the cleaning step
fit = obj.clean(fit);

end


%% LOCAL FUNCTIONS

function hrf = makeGammaHRF(x, stimDeltaT)

% Construct an HRF
gamma1 = x(1);
gamma2 = x(2);
undershootGain = x(3);
duration = 28; % Fixed value in seconds to match FLOBS vector lengths

% Define a timebase at the data resolution
timebase = 0:stimDeltaT:duration;

% Create the double gamma function
g1 = gampdf(timebase,gamma1, 1);
g1 = g1./ max(g1);
g2 = gampdf(timebase, gamma2, 1);
g2 = (g2/ max(g2)) * undershootGain;
hrf = g1 - g2;

% Set to zero at onset
hrf = hrf - hrf(1);

% Normalize the kernel to have unit area.
hrf = hrf/sum(abs(hrf));

% Make the hrf a column vector
hrf = hrf';

end


function hrf = makeFlobsHRF(x, flobsbasis)

% Create the HRF
hrf = flobsbasis*x';

% Normalize the kernel to have unit area
hrf = hrf/sum(abs(hrf));

end

