function x0 = initial(obj)
% Returns initial guess for the model parameters
%
% Syntax:
%   x0 = obj.initial()
%
% Description:
%   Initial values for the prf_timeShift model. Rationale is as follows:
% - gain
% - exponent
% - cardinal multiplier
% - time-constant
%
% Inputs:
%   none
%
% Optional key/value pairs:
%   none
%
% Outputs:
%   x0                    - 1xnParams vector.
%


% Obj variables
nParams = obj.nParams;
nFilterBins = obj.nFilterBins;
nFixedParamsAdapt = obj.nFixedParamsAdapt;
mu = obj.mu;

% Assign the x0 variable
x0 = zeros(1,nParams);

% Assemble X0
x0(1) = 0; % gain
x0(2) = 1; % exponent
x0(3) = 1; % tau=1 (in seconds) of the exponential integrator of heading direction

% zero initial gain for direction filter bank
x0(nFixedParamsAdapt+1:nFixedParamsAdapt+nFilterBins) = 0;

% HRF params
x0(nParams-2:nParams) = mu; % FLOBS eigen1, 2, 3

end

