 % crossValDemo.m
% tbUseProject('headingModel')
% This script demonstrates working with model objects for the purpose of
% performing cross-validation analyses.


% The mat file we need is located in the "data" directory. Construct a path
% to this file, starting with the location of this script.
parcels = {'bMask', 'EVC', 'OPA', 'PPA', 'RSC', 'PHG', 'ERC', 'Hipp', 'Thal'...
    'VTC'};
for sub_n = 8:8
    if sub_n <10
        sub = ['sub-0' num2str(sub_n)];
    else
        sub = ['sub-' num2str(sub_n)];
    end
    % sub='sub-08';
    for p=5:5
        parcel=parcels{p};
%         parcel='PHG';
        fileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),['data/' sub], ...
        [sub '_city1A_stimulus_data_' parcel '-city1AB.mat']);

        % Load the stimulus and data variables
        load(fileName,'stimulus','data')

        % The TR of the experiment, in seconds
        tr=2;

        % Define the set of vertices
        vxs = 1:size(data{1},1);

        % Let's just work with one vertex for the demo
        % ii = 1;
        % vxs = vxs(ii);

        % Define modelOpts
        modelOpts = {'nFilterBins',45};

        % Call the forwardModel
        results = forwardModel(data,stimulus,tr,'modelClass','heading','vxs',vxs,'modelOpts',modelOpts);

        %% VALIDATION
        % Note that the steps below will over-write the "stimulus" and "data"
        % variables. We may wish to make this more robust later.

        % Load the validation dataset
        fileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),['data/' sub],...
            [sub '_city1B_stimulus_data_' parcel '-city1AB.mat']);
        load(fileName,'stimulus','data')

        % Re-create the model, now with the new stimulus and data
        model = heading(data,stimulus,tr,modelOpts{:});

        % Prepare the validation data into a signal. You would want to loop over
        % vertices to test this across the entire dataset.
        dataPrep = catcell(2,model.prep(data));
        % ii = 1; % Just do the first vertex
        cv_R2=zeros(length(vxs), 1);
        for ii=1:length(vxs)
            signal = model.clean(dataPrep(vxs(ii),:)');

            % Obtain the parameter values of the previous fit for this vertex
            x = results.params(vxs(ii),:);

            % Obtain the R2 metric for the signal, given the params (x) we have
            % previously calculated
            valR2 = model.metric(signal,x);
            cv_R2(ii, 1)=valR2;
        end

        % save R2
        outFileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),'results',...
            [sub '_city1A_forward-heading-45_city1B_' parcel '-city1AB.mat']);
        save(outFileName,'cv_R2','results');
    end
end
