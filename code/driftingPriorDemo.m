

%% Fixed variables of the simulation
%   - tr in seconds
%   - The preferred direction to be modeled, in radians
%   - The number of TRs to extend the first heading value prior to the
%     start ofthe scan. This is necessary to avoid a ramp up of BOLD fMRI
%     response at the start of the experiment. This could also be
%     accomplished by discarding the first n TRs of the data in the
%     modeling effort.
tr = 2;
preferredDirection = pi;
nFixedParams = 3; % corresponding to the adaptation gain, epsilon, and muu

% The temporal sampling interval (in seconds) of the heading direction
% vector
deltaTStim = 660/3384;

% The number of seconds to pad the start of the stimulus sequence, to
% provide some "warm up" for the model
padTimeSecs = 8*2;
nDTsToPad = round(padTimeSecs*(1/deltaTStim));

%% Create a stimulus
% We load a set of heading direction vectors from an example subject, and
% either use these values, or base a simulation upon them.
fileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),'data/', ...
    'sub-08_city1A_stimulus.mat');
load(fileName,'stimulus')
nDTs = size(stimulus{1},2);
nTRs = 330; % TRs per acquisition

% Extend the stimulus back in time to provide some "warm up" for the HRF;
% create a stimTime array
stimTime={};
for ii=1:length(stimulus)

    % Grab an acquisition from the cell array
    thisStim = stimulus{ii};

    % Create a stimTime vector
    stimTime{ii} = -nDTsToPad*deltaTStim:deltaTStim:deltaTStim*(length(thisStim)-1);

    % Pad this stim
    thisStim = [repmat(thisStim(1),1,nDTsToPad), thisStim];
    stimulus{ii} = thisStim;

end

% Do the demo for one acquisition
stimulus = stimulus{1};

% Create a figure to show the results
figure

% These are the mu values that we will show
muVals = [0,0.15:0.2:0.95];

for mm = 1:length(muVals)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% THIS IS THE SECTION THAT WOULD GO INTO THE FORWARD MODEL

    % Define an empty heading change vector
    headingChange = zeros(size(stimulus));

    headingPrior = zeros(size(stimulus));

    % Set the initial value of the heading prior. If mu==0, then the prior
    % never moves, so put it in the middle of the stimulus range.
    if mm==1
        headingPrior(1)=0;
    else
        headingPrior(1)=stimulus(1);
    end

    % Loop through the subsequent events and update the prior
    for n=2:length(stimulus)
        angChange=wrapTo2Pi(stimulus(n)-stimulus(n-1));
        headingPrior(n)=wrapTo2Pi(headingPrior(n-1)+(1-muVals(mm))*angChange);
    end

    adaptEffect = wrapTo2Pi(stimulus-headingPrior);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Plot the results
    subplot(length(muVals),2,(mm-1)*2+1)
    plot(headingPrior,'-k'); hold on
    title(sprintf('heading Prior, mu = %2.2f',muVals(mm)))
    ylim([0 2*pi]);
    box off

    subplot(length(muVals),2,(mm-1)*2+2)
    plot(adaptEffect,'-r'); hold on
    title(sprintf('adaptation effect, mu = %2.2f',muVals(mm)))
    ylim([0 2*pi]);
    box off

end