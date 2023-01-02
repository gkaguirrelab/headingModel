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

% Obj variables
stimAcqGroups = obj.stimAcqGroups;
stimTime = obj.stimTime;
nParams = obj.nParams;
nFixedParamsAdapt = obj.nFixedParamsAdapt;
filterResponse = obj.filterResponse;

% Break the parameters of x into named variables for code transparency
adaptGain = x(1);   % Gain of the adaptation effect
epsilon = x(2);     % Non-linear exponent of the neural signal
tau = x(3);         % Time constant of the temporal integration


%% Build a model of absolute heading direction
% Use the stored filter bank from the object properties.
% Apply the filter weights as a single matrix multiplication.
filterWeights=x(nFixedParamsAdapt+1:end-3);
neuralSignal = sum(filterResponse.*filterWeights, 2);


%% Load or create the headingChange
% Check if the tau parameter has changed, and load or generate the
% headingChange as appropriate
if tau == obj.tauLast

    headingChange = obj.headingChangeLast;

else

    % Obtain the temporal support for one the stimulus for one acquisition
    stimTimeSingleAcq = stimTime(stimAcqGroups==1);

    % Create an exponential kernel under the contol of tau (in units of
    % seconds)
    exponentialKernel = exp(-1/tau*stimTimeSingleAcq);

    % Normalize the kernel to have unit area
    exponentialKernel = exponentialKernel/sum(exponentialKernel);

    % Convolve the unwrapped stimulus by the exponential kernel
    headingPrior = conv2run(obj.stimulusUnwrap,exponentialKernel,stimAcqGroups);

    % The heading change is simply the (unwrapped) angular distance between the
    % current heading and the integrated prior heading
    headingChange = abs(obj.stimulusUnwrap - headingPrior);

    % Store the tau and headingPrior values in the object
    obj.headingChangeLast = headingChange;
    obj.tauLast = tau;

end

% Apply an exponential parameter to produce a compressive or expansive
% non-linear mapping between heading change and neural response (a
% parameter value of 1 provides a linear mapping).
headingChange = headingChange.^epsilon;

% Scale the stimulus matrix by the gain parameter
neuralSignal = neuralSignal + headingChange*adaptGain;

% Create the HRF 
if all(x(nParams-2:nParams)==0) % We want a delta function
    hrf = zeros(size(obj.flobsbasis,1),1);
    hrf(1)=1;
else
    hrf = makeFlobsHRF(x(nParams-2:nParams), obj.flobsbasis);
end

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


function hrf = makeFlobsHRF(x, flobsbasis)

% Create the HRF
hrf = flobsbasis*x';

% Normalize the kernel to have unit area
hrf = hrf/sum(abs(hrf));

end

