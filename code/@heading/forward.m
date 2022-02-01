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
stimulus = obj.stimulus;
stimAcqGroups = obj.stimAcqGroups;
stimTime = obj.stimTime;
nParams = obj.nParams;
dataTime = obj.dataTime;

%% IMPLEMENT THE NON-LINEAR MODEL OF HEADING RESPONSE HERE
% In the code here, you would take the stimulus vector, and perform
% manipulations upon it to yield the time-series prediction. For example,
% you might
% 1) Take the first-derivative of the continuous measure of heading
% direction, and then ciruclarize it so that the vector reflects the
% absolute magnitude of heading change in any given TR
% 2) Scale this vector by the first parameter (to implement a "gain"
% effect)
% 3) Convolve the vector by an exponential decay kernel, with the
% time-constant of the exponential decay under the control of the second
% parameter of the model.
% 4) Convolve the resultant vector by the HRF, with the shape of the hRF
% under the control of the final 3 parameters of the model.

% Break the parameters of x into named variables for code transparency
gain = x(1);
epsilon = x(2);
beta = x(3);
tau = x(4);

% Define an empty heading change vector
headingChange = zeros(size(stimulus),class(stimulus));

% Loop over acquisitions and obtain the absolute, circular value of the
% first derivative of the heeading vector
for run = 1:max(stimAcqGroups)
    temp = stimulus(stimAcqGroups==run);
    headingChange(stimAcqGroups==run,:) = [0;abs(angdiff(temp))];
end

% Apply an exponential parameter to produce a compressive or expansive
% non-linear mapping between heading change and neural response (a
% parameter value of 1 provides a linear mapping).
headingChange = headingChange.^epsilon;

%% COULD IMPLEMENT HERE THE USE OF BETA TO ENHANCE THE EFFECT OF THE
% cadinal meridians upon the neurall response

% Obtain the temporal support for one acquisition
dataTimeSingleAcq = dataTime(stimAcqGroups==1);

% Create an exponential kernel under the contol of tau (in units of
% seconds)
exponentialIRF = exp(-1/tau*dataTimeSingleAcq);

% Normalize the kernel to have unit area
exponentialIRF = exponentialIRF/sum(abs(exponentialIRF));

% Apply the exponential kernel to the heading time series
headingChange = conv2run(headingChange,exponentialIRF,stimAcqGroups);

% Scale the stimulus matrix by the gain parameter
neuralSignal = headingChange*gain;

% Create the HRF
switch obj.hrfType
    case 'flobs'
        hrf = makeFlobsHRF(x(nParams-2:nParams), obj.flobsbasis);
    case 'gamma'
        hrf = makeGammaHRF(x(nParams-2:nParams), obj.stimDeltaT);
    otherwise
        error('Not a recognized HRF type')
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


