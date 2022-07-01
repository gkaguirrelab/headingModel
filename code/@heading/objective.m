function fVal = objective(obj, signal, x)
% Evaluates the match between a signal and a model fit
%
% Syntax:
%   fVal = obj.objective(signal, x)
%
% Description:
%   Given a time series signal and the parameters of the forward model,
%   returns the objective function to be minimized. In some applications,
%   this could just be the negative of obj.metric.
%
% Inputs:
%   signal                - 1 x time vector. The data to be fit.
%   x                     - 1 x nParams vector of parameter values.
%
% Optional key/value pairs:
%   none
%
% Outputs:
%   fVal                  - Scalar.
%

% Implement an L2 norm
fVal = std(signal - obj.forward(x));

% Add a lasso regression penalty based upon the bin weights
lassoRegularization = obj.lassoRegularization;
nFilterBins = obj.nFilterBins;
nFixedParams = obj.nFixedParams;
penalty = lassoRegularization * sum(abs(x(nFixedParams+1:nFixedParams+nFilterBins)))/nFilterBins;
fVal = fVal + penalty;

end

