% Set some fixed parameters, and see if we can recover them, including
% an interpolated preferred heading direction.
%% set up
noise=0.5;
gain=0.25;
epsilon=0.8;
tau_range=0.5:0.5:10;
nBins = 45;
nSim=25;
gain_re=zeros(length(tau_range), nSim);
epsilon_re=zeros(length(tau_range), nSim);
tau_re=zeros(length(tau_range), nSim);
for t = 1:length(tau_range)
    tau=tau_range(t);
    for n=1:nSim
        fixedParams = [gain epsilon tau];
        useRealHeading = true; hrfSearch = false; adaptSearch = true;
        binWts=0; lassoRegularization=0.05;
        [x0, x1] = simEngine(noise,binWts,fixedParams,lassoRegularization,...
            nBins,nBins,useRealHeading,hrfSearch,adaptSearch,1);
        gain_re(t,n)=x1(1);
        epsilon_re(t,n)=x1(2);
        tau_re(t,n)=x1(3);
    end
end

figure
subplot(3,1,1);
plot(tau_range, mean(y1,2));
ylabel('Recovere Gain');
% set(gca,'YTick',0:pi/2:2*pi)
% set(gca,'YTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'})
% xlim([min(headingTime) max(headingTime)]);

subplot(3,1,2);
plot(tau_range, mean(y2,2));
ylabel('Recovere Epsilon');
% title(sprintf('simulated neural response downsampled to TRs (bins = %d)',nSimBins));
% xlim([min(dataTime) max(dataTime)]);

subplot(3,1,3);
plot(tau_range, mean(y3,2));
xlabel('Tau');
ylabel('Recovere Tau');
% title(sprintf('simulated and fitted BOLD response (bins = %d)',nSimBins));
% xlim([min(dataTime) max(dataTime)]);

y1=gain_re-gain;
y2=epsilon_re-epsilon;
y3=tau_re-tau_range';
figure();
plot(tau_range, mean(y1,2));
hold on
plot(tau_range, mean(y2,2));
hold on
plot(tau_range, mean(y3,2));
figure
for i=1:length(tau_range)
    scatter(tau_re(i, :))
    hold on
end