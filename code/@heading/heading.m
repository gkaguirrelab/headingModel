classdef heading < handle
    
    properties (Constant)
        
        % Properties of the search stages.
        nStages = 1;
        
        % A description of the model
        description = ...
            ['A model of neural response to change in heading.\n'];
    end
    
    % Private properties
    properties (GetAccess=private)
        
        % The FLOBS eigenvectors
        flobsbasis
        
        % The multivariate normal means of the 3 eigenvectors
        mu
        
        % The multivariate normal covariance matrix
        C
        
        % The projection matrix used to regress our nuisance effects
        T
        
    end
    
    % Calling function can see, but not modify
    properties (SetAccess=private)
        
        % The number of bins in the heading direction filter bank
        nFilterBins

        % The number of parameters dedicated to the adaptation model
        nFixedParamsAdapt

        % The number of parameters in the model
        nParams
        
        % Properties of the search stages.
        floatSet
        fixSet
        
        % A vector of the length totalTRs x 1 that has an index value to
        % indicate which acquisition (1, 2, 3 ...) a data time
        % sample is from.
        dataAcqGroups
        
        % A vector of length totalTRs x 1 and in units of seconds that
        % defines the temporal support of the data relative to the time of
        % onset of the first TR of each acquisition.
        dataTime
        
        % TR of the data in seconds
        dataDeltaT
        
        % The stimulus vector, concatenated across acquisitions and
        % squished across x y. Thus it will have the dimensions:
        %	[totalST x*y]
        stimulus
        
        % A set of labels for the stimuli, used to label maps and the
        % result fields
        stimLabels
        
        % A vector of length totalST x 1 that has an index value to
        % indicate which acquisition (1, 2, 3 ...) a stimulus time
        % sample is from.
        stimAcqGroups
        
        % A vector of length totalST x 1 and in units of seconds that
        % defines the temporal support of stimulus relative to the time of
        % onset of the first TR of each acquisition. If set to empty, the
        % stimTime is assumed to be equal to the dataTime.
        stimTime
        
        % The temporal resolution of the stimuli in seconds.
        stimDeltaT
        
        % A cell array that contains things that the model might want
        payload
        
        % A time x 1 vector that defines the HRF convolution kernel
        hrf
        
        % The lasso regression penalty for the bin weights
        lassoRegularization

        % The unwrapped stimulus
        stimulusUnwrap
        
        % The filter response, using circular gaussian function
        filterResponse

        % We store the tauLast and headingChangeLast to avoid
        % re-computing these if tau has not changed
        tauLast
        headingChangeLast

    end
    
    % These may be modified after object creation
    properties (SetAccess=public)
        
        % The number of low frequencies to be removed from each acquisition
        polyDeg
        
        % Typical amplitude of the BOLD fMRI response in the data
        typicalGain
        
        % The lower and upper bounds for the model
        lb
        ub
        
        % A vector, equal in length to the number of parameters, that
        % specifies the smallest step size that fmincon will take for each
        % parameter. This threshold is also used to determine if, in a call
        % to obj.forward, the gaussvector needs to be re-calculated, or if
        % the prior value can be used.
        FiniteDifferenceStepSize
        
        % Verbosity
        verbose
        
    end
    
    methods
        
        % Constructor
        function obj = heading(data,stimulus,tr,varargin)
            
            % instantiate input parser
            p = inputParser; p.KeepUnmatched = false;
            
            % Required
            p.addRequired('data',@iscell);
            p.addRequired('stimulus',@iscell);
            p.addRequired('tr',@isscalar);
            
            p.addParameter('stimTime',{},@iscell);
            p.addParameter('stimLabels',{},@iscell);
            p.addParameter('payload',{},@iscell);
            p.addParameter('polyDeg',[],@isnumeric);
            p.addParameter('typicalGain',300,@isscalar);
            p.addParameter('nFilterBins',8,@isscalar);
            p.addParameter('adaptSearch',true,@islogical);            
            p.addParameter('hrfSearch',true,@islogical);            
            p.addParameter('lassoRegularization',0.05,@isscalar);           
            p.addParameter('verbose',true,@islogical);
            
            % parse
            p.parse(data, stimulus, tr, varargin{:})
            
            % Store the lassoRegularization value
            obj.lassoRegularization = p.Results.lassoRegularization;

            % Create the dataTime and dataAcqGroups variables
            % Concatenate and store in the object.
            for ii=1:length(data)
                dataAcqGroups{ii} = ii*ones(size(data{ii},2),1);
                dataTime{ii} = (0:tr:tr*(size(data{ii},2)-1))';
            end
            obj.dataAcqGroups = catcell(1,dataAcqGroups);
            obj.dataTime = catcell(1,dataTime);
            obj.dataDeltaT = tr;
            clear data
            
            % Each row in the stimulus is a different stim type that will
            % be fit with its own gain parameter. Record how many there are
            nStimTypes = size(stimulus{1},1);
            
            % These are the parameters, corresponding to:
            % - gain
            % - exponent
            % - tau
            % - a variable number of parameters for an absolute heading direction model
            % - 3 parameters of the FLOBS HRF
            obj.nFixedParamsAdapt = 3;
            obj.nFilterBins = p.Results.nFilterBins;
            nHRFParams = 3;
            obj.nParams = obj.nFixedParamsAdapt + p.Results.nFilterBins + nHRFParams;

            % Define the stimLabels
            if ~isempty(p.Results.stimLabels)
                stimLabels = p.Results.stimLabels;
                if length(stimLabels) ~= nStimTypes
                    error('forwardModelObj:badStimLabels','the stimLabels value must be a cell array equal to the number of stimulus types.');
                end
            else
                stimLabels = cell(1,nStimTypes);
                for pp = 1:nStimTypes
                    stimLabels{pp} = sprintf('beta%02d',pp);
                end
            end
            obj.stimLabels = stimLabels;
            
            % Define the fix and float param sets. We will always search
            % over the nFixedParamsOther and the nFilterBins
            obj.fixSet = {1:obj.nParams-3, 1:obj.nParams-3};
            obj.floatSet = {obj.nParams-2:obj.nParams, obj.nParams-2:obj.nParams};

            floatSet = [obj.nFixedParamsAdapt+1:obj.nFixedParamsAdapt+obj.nFilterBins];
            fixSet = [];

            if p.Results.hrfSearch
                floatSet = [floatSet obj.nParams-nHRFParams+1:obj.nParams];
            else
                fixSet = [fixSet obj.nParams-nHRFParams+1:obj.nParams];
            end

            if p.Results.adaptSearch
                floatSet = [floatSet 1:obj.nFixedParamsAdapt];
            else
                fixSet = [fixSet 1:obj.nFixedParamsAdapt];
            end

            obj.floatSet = {floatSet};
            obj.fixSet = {fixSet};

            % Create the stimulus and stimAcqGroups variables. Concatenate
            % the cells and store in the object.
            for ii=1:length(stimulus)
                % Transpose the stimulus matrix within cells
                stimulus{ii} = stimulus{ii}';
                stimAcqGroups{ii} = ii*ones(size(stimulus{ii},1),1);
            end            
            obj.stimulus = catcell(1,stimulus);
            obj.stimAcqGroups = catcell(1,stimAcqGroups);
            
            % Create an unwrapped version of the stimulus to be used later.
            % Do so while respecting acquisition boundaries
            stimulusUnwrap = nan(size(obj.stimulus));
            for ii=1:length(stimulus)
                stimulusUnwrap(obj.stimAcqGroups==ii) = unwrap(stimulus{ii});
            end
            obj.stimulusUnwrap = stimulusUnwrap;

            % Construct and / or check stimTime
            if isempty(p.Results.stimTime)
                % If stimTime is empty, check to make sure that the length
                % of the data and stimulus matrices in the time domain
                % match
                if length(obj.stimAcqGroups) ~= length(obj.dataAcqGroups)
                    error('forwardModelObj:timeMismatch','The stimuli and data have mismatched temporal support and no stimTime has been passed.');
                end
                % Set stimTime to empty
                obj.stimTime = [];
                % The temporal resolution of the stimuli is the same as the
                % temporal resolution of the data
                obj.stimDeltaT = tr;
            else
                % We have a stimTime variable.
                stimTime = p.Results.stimTime;
                % Make sure that all of the stimTime vectors are regularly
                % sampled (within 3 decimal precision)
                regularityCheck = cellfun(@(x) length(unique(round(diff(x),3))),stimTime);
                if any(regularityCheck ~= 1)
                    error('forwardModelObj:timeMismatch','One or more stimTime vectors are not regularly sampled');
                end
                % Make sure that the deltaT of the stimTime vectors all
                % match, and store this value
                deltaTs = cellfun(@(x) x(2)-x(1),stimTime);
                if length(unique(deltaTs)) ~= 1
                    error('forwardModelObj:timeMismatch','The stimTime vectors do not have the same temporal resolution');
                end
                obj.stimDeltaT = deltaTs(1);
                % Concatenate and store the stimTime vector
                for ii=1:length(stimTime)
                    % Transpose stimTime within cells
                    stimTime{ii} = stimTime{ii}';
                end
                obj.stimTime = catcell(1,stimTime);
                % Check to make sure that the length of the stimTime vector
                % matches the length of the stimAcqGroups
                if length(obj.stimTime) ~= length(obj.stimAcqGroups)
                    error('forwardModelObj:timeMismatch','The stimTime vectors are not equal in length to the stimuli');
                end
            end

            % Create and store the bank of filters that will be used to
            % model absolute heading direction
            binSeparation = (2*pi/obj.nFilterBins);
            binCenters = 0:binSeparation:(2*pi)-binSeparation;
            sigma = binSeparation/(2*sqrt(2*log(2)));
            kappa = 1/sigma^2;
            obj.filterResponse=circ_vmpdf(obj.stimulus,binCenters,kappa);

            % Initialize the tau and exponential kernel storage
            obj.tauLast = nan;
            obj.headingChangeLast = nan;

            % Done with these big variables
            clear data stimulus stimTime acqGroups
            
            % Distribute other params to obj properties
            obj.payload = p.Results.payload;
            obj.polyDeg = p.Results.polyDeg;
            obj.typicalGain = p.Results.typicalGain;
            obj.verbose = p.Results.verbose;
            
            % Create and cache the flobs basis
            obj.genflobsbasis;
            
            % Set the bounds and minParamDelta
            obj.setbounds;
            
            % Create and cache the projection matrix
            obj.genprojection;
            
            % Call the forward model to create and store an initial hrf
            obj.forward(obj.initial);
            
            
        end
        
        % Required methds -- The forwardModel function expects these
        rawData = prep(obj,rawData)
        x0 = initial(obj)
        signal = clean(obj, signal)
        [c, ceq] = nonlcon(obj, x)
        fVal = objective(obj, signal, x)
        [fit, hrf] = forward(obj, x)
        x0 = update(obj,x,x0,floatSet,signal)
        metric = metric(obj, signal, x)
        seeds = seeds(obj, data, vxs)
        results = results(obj, params, metric)
        results = plot(obj, data, results)
        
        % Internal methods
        genflobsbasis(obj);
        setbounds(obj)
        genprojection(obj)
        
        % Set methods
        function set.polyDeg(obj, value)
            obj.polyDeg = value;
            obj.genprojection;
        end
        
    end
end