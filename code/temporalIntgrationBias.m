%% calculates, for each point in time, the “prior” of current heading 
%% direction (prior in this sense meaning an integrated average over some 
%% time history).
muRange = 0:0.1:0.9;
nMus=length(muRange);
%% Obtain the vector of heading directions (not heading change; the vector 
% of raw heading directions)
fileName = fullfile(fileparts(fileparts(mfilename('fullpath'))),'data/', ...
    'sub-08_city1A_stimulus.mat');
load(fileName,'stimulus')
%% Loop over time points. At each time point, calculated an updated value 
% of r. To do so, you will take the previous value of r, and add to it 
% the weighted difference between that prior value and the current heading 
% direction. This difference needs to be the circular difference
% headingPrior={};
run=1;
headingRun=stimulus{run};
% headingRun=180*headingRun/pi;
% nTP=length(headingRun);
nTP=1200;
figure
subplot(1+nMus,1,1);
plot(1:nTP,headingRun(1:nTP),'-ok');
% ylim([0 360]);
ylim([0 2*pi]);
title('Heading');
hold on
for m=1:nMus
    mu=muRange(m);
% for run=1:length(stimulus)
    % Set the initial state of the prior (rzero) to be the same as the heading 
    % direction for the first time point
    r0=headingRun(1);
    
    headingPriorRun = zeros(size(headingRun));
    headingPriorRun(1)=r0;
    for n=2:nTP
        currentHeading=headingRun(n);
        angChange=angdiff(r0,currentHeading);
%         angChange=pi+angChange;
%         if angChange<0
%             angChange=2*pi+angChange;
%         end
%         angChange=mod(angChange,2*pi);
        r=r0+(1-mu)*angChange;
%         if r>2*pi
%             r=mod(r,2*pi);
%         end
        r0=r;
        headingPriorRun(n)=r;
    end
    subplot(1+nMus,1,m+1);
    headingPriorRun=wrapTo2Pi(headingPriorRun);
    plot(1:nTP,headingPriorRun(1:nTP),'-ok');
    ylim([0 2*pi]);
    title(['Heading Prior (mu: ' num2str(mu) ')']);
    
%     headingPrior{run}=headingPriorRun;
% end
end
%% Create a plot that shows the heading directions, followed by a plot 
% that shows the calculated prior for different values of the mu parameter.
