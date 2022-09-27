function comp_ln_nl(voxel_index)
%% which 
    tbUseProject('headingModel')
    clear;
    clc;
    %% Load the weights from different models
    % load nonlinear model with false adapt flag set
    load('./results/sub-08_city1A_head-45_city1B_bMask_nl-lasso-02_8TR-pad_adapt-0.mat', 'model_params');
    nl_adpt_0_wts = model_params;
%     nl_adpt_0_r2 = cv_R2;
    % load nonlinear model with true adapt flag
    load('./results/sub-08_city1A_head-45_city1B_bMask_nl-lasso-02_8TR-pad.mat', 'model_params')
    nl_wts=model_params;
%     nl_r2 = cv_R2;
    % load linear model
    load('./results/sub-08_city1A_head-45_city1B_bMask_ln.mat', 'train_weights', 'test_DM')
    ln_wts=train_weights;
%     ln_r2=test_R2;
    %% VALIDATION
    % Note that the steps below will over-write the "stimulus" and "data"
    % variables. We may wish to make this more robust later.
    sub='sub-08';
    parcel='bMask';
    test_dst = 'B';
    nTRsToPad = 8;
    tr=2;
    % Define modelOpts
    nFitBins = 45;
    hrfSearch = true;
    adaptSearch = false;
    lassoRegularization=0;
    modelOpts = {'nFilterBins',nFitBins,'hrfSearch',hrfSearch,'adaptSearch',adaptSearch,'typicalGain',1, ...
    'lassoRegularization',lassoRegularization};
    % Load the validation dataset
    fileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),['data/' sub],...
        [sub '_city1' test_dst '_stimulus_data_' parcel '.mat']);
    load(fileName,'stimulus','data')

    % Extend the stimulus back in time to provide some "warm up" for the HRF;
    % create a stimTime array
    for ii=1:length(stimulus)
        % Grab an acquisition from the cell array
        thisStim = stimulus{ii};
        % Create a stimTime vector
        test_stimTime{ii} = -nTRsToPad*tr:tr:tr*(length(thisStim)-1);
        % Pad this stim
        thisStim = [repmat(thisStim(1),1,nTRsToPad), thisStim];
        stimulus{ii} = thisStim;
    end
    % Re-create the model, now with the new stimulus and data
    model = heading(data,stimulus,tr,'stimTime',test_stimTime,modelOpts{:});
    % Prepare the validation data into a signal. You would want to loop over
    % vertices to test this across the entire dataset.
    dataPrep = catcell(2,model.prep(data));
    %% Define the set of vertices
    vxs = 1:size(data{1},1);
    % Just do this vertex or vertex which showed the best r2 in ln model
    if voxel_index == 0
        [max_r2, ii]=max(ln_r2)
    elseif which_vox == -1
        [max_r2, ii]=max(nl_r2)
    else
        ii = voxel_index;
    end
    signal = model.clean(dataPrep(vxs(ii),:)');
    % Obtain the parameter values of the previous fit for this vertex
    nl_adpt_0_wts_this_vox=nl_adpt_0_wts(vxs(ii),:);
    nl_wts_this_vox=nl_wts(vxs(ii),:);
    ln_wts_this_vox=ln_wts(vxs(ii),:);
    nl_ln_combined_this_vox=[nl_adpt_0_wts_this_vox(1:3), ln_wts_this_vox, ...
        nl_adpt_0_wts_this_vox(end-2:end)];
    %% calculate the predicted time course yhat
    yhat_nl_adpt_0 = model.forward(nl_adpt_0_wts_this_vox);
    yhat_nl = model.forward(nl_wts_this_vox);
    yhat_nl_ln_combined = model.forward(nl_ln_combined_this_vox);
    % Obtain the R2 metric for the signal, given the params (x) we have
    % previously calculated
    valR2_nl_adpt_0 = model.metric(signal,nl_adpt_0_wts_this_vox);
    valR2_nl = model.metric(signal,nl_wts_this_vox);
    valR2_nl_ln_combined = model.metric(signal,nl_ln_combined_this_vox);
    %% calculate the linear model yhat
    test_data=data;
    test_tCourse = cat(2,test_data{:})';
    yhat_ln=sum(repmat(ln_wts_this_vox,numel(test_tCourse(:,ii)),1).*[ones(1,numel(test_tCourse(:,ii))); test_DM]',2);
    valR2_ln=(corr(test_tCourse(:,ii),test_yhat, 'type', 'Pearson'))^2;

    %% Make plots
    x=1:length(yhat_ln)/9;
    figure
    subplot(4,1,1);
    plot(x,yhat_ln(x), 'b');
    hold on
    plot(x,test_tCourse(x,ii),'-k');
    ylabel('Neural response');
    title(sprintf('Linear ridge r2: %s', num2str(valR2_ln)));
    legend('predicted', 'real BOLD');
    xlim([1 330]);
    
    subplot(4,1,2);
    plot(x,yhat_nl_ln_combined(x),'-r');
    hold on
    plot(x,signal(x),'-k');
    ylabel('Neural response');
    title(sprintf('Linear Nonlinear parameters combined r2: %s', num2str(valR2_nl_ln_combined)));
    legend('predicted', 'real BOLD');
    xlim([1 330]);

    subplot(4,1,3);
    plot(x,yhat_ln(x),'b');
    hold on
    plot(x,yhat_nl_ln_combined(x), 'r');    
    ylabel('Neural response');
    title('Predicted BOLD: linear vs. linear + nonlinear');
    xlim([1 330]);
    legend('linear', 'linear+nonlinear');
    
    subplot(4,1,4);
    plot(x,yhat_nl(x),'g'); 
    hold on
    plot(x,yhat_nl_ln_combined(x), 'r');
    ylabel('Neural response');
    title('Predicted BOLD: nonlinear vs. linear + nonlinear');
    xlim([1 330]);
    legend('nonlinear', 'linear+nonlinear');
    xlabel('TRs (One Run)');
%     ylabel('weights');
% xlabel('time [seconds]');
    % set(gca,'YTick',0:pi/2:2*pi)
    % set(gca,'YTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'})
    % xlim([min(dataTime) max(dataTime)]);
end