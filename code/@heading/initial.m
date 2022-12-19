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
% nFixedParamsOther = obj.nFixedParamsOther;

% Assign the x0 variable
x0 = zeros(1,nParams);

% Assemble X0
x0(1) = 0; % gain
x0(2) = 1; % exponent
x0(3) = 0; % tau=1 (in seconds) of the exponential integrator of heading direction
           % muu=0 no memory at all
x0(nFixedParamsAdapt+1:nFixedParamsAdapt+nFilterBins) = 0;       % zero initial gain for direction filter bank

switch obj.hrfType
    case 'flobs'
        x0(nParams-2:nParams) = [0.86, 0.09, 0.01]; % FLOBS eigen1, 2, 3
    case 'gamma'
        x0(nParams-2:nParams) = [6, 10, 0.1]; % Gamma params
    otherwise
        error('Not a valid hrfType')
end

end

