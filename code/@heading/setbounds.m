function setbounds(obj)
% Sets the bounds on the model parameters
%
% Syntax:
%   obj.setbounds()
%
% Description:
%   Bounds for the prf_timeShift model. Rationale is as follows:
% - gain
% - exponent
% - cardinal multiplier
% - time-constant
%
%   These are specified as 1 x nParams vectors.
%
% Inputs:
%   none
%
% Optional key/value pairs:
%   none
%
% Outputs:
%   none
%

% Obj variables
nParams = obj.nParams;

% Define outputs
lb = nan(1,nParams);
ub = nan(1,nParams);

% The gain parameter is unbounded
lb(1) = -Inf;             % gain
ub(1) = Inf;              % gain

% Raise the heading change to an exponent to support compressive /
% expansive non-linearities
lb(2) = 0.1;             % gain
ub(2) = 3;              % gain

% Enhance the value of direction changes close to the cardinal meridians
lb(3) = 0.5;             % gain
ub(3) = 3;              % gain

% The time-constant parameter is bounded by zero at the low end, and by 2
% seconds at the high end. We set 2 as the upper bound to avoid colliding
% with the HRF model
lb(4) = 0.01;
ub(4) = 2;

% The HRF shape parameters vary by model type
switch obj.hrfType
    case 'flobs'
        
        % Object properties associated with the FLOBS eigenvectors
        mu = obj.mu;
        C = obj.C;
        
        % Set bounds at +-15SDs of the norm distributions of the FLOBS
        % parameters
        sd15 = 15*diag(C)';
        
        lb(nParams-2:nParams) = mu-sd15;	% FLOBS eigen1, 2, 3
        ub(nParams-2:nParams) = mu+sd15;	% FLOBS eigen1, 2, 3

    case 'gamma'
        lb(nParams-2:nParams) = [2 6 0];	% Gamma1,2, and undershoot gain
        ub(nParams-2:nParams) = [8 12 2];	% Gamma1,2, and undershoot gain

    otherwise
        error('Not a valid hrfType')
end

% Store the bounds in the object
obj.lb = lb;
obj.ub = ub;

% Store the FiniteDifferenceStepSize for the model. See here for more
% details:
%   https://www.mathworks.com/help/optim/ug/optimization-options-reference.html
FiniteDifferenceStepSize = nan(1,nParams);
FiniteDifferenceStepSize(1,:) = sqrt(eps);
obj.FiniteDifferenceStepSize = FiniteDifferenceStepSize;

end

