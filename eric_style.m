N=25;
theta = zeros(1,2*N+1)+0.001*2*pi*randn(1,2*N+1);
locations=-N:N;
Tfinal = 18000;

% % base period varies gradually and symmetrically (linearly or quadratically) across the jaw
E0 = 30;
E1 = E0/14;
%E = E0 - E1*(locations/N).^2;
E = E0 + E1*abs(locations/N);
% E = E0 - 2.7*abs(locations/N);
% E=E0;
m = 1;          % sharpness of signal influence as a function of phases
epsilon = 0.1; % strength of coupling
dt = 0.1;
Nsteps=Tfinal/dt;
theta = repmat(theta,Nsteps+1,1); %preallocate solution for better performance


for k = 1:Nsteps
    theta(k+1,:) = theta(k,:) + dt* (  2*pi./E - ...
        epsilon * ( [0 ( (1+cos(theta(k,1:end-1)) )/2 ).^m ]+ ...
        [ ( ( 1+cos(theta(k,2:end)))/2).^m 0] ) .* ((1-sin(theta(k,:)))/2).^m );
end
[i,j] = find(diff(mod(theta,2*pi))<0);
EruptionTimes = [locations(j)',i*dt];
EruptionTimes = sortrows(EruptionTimes,2);

%%
figure(4);
imagesc(mod(theta,2*pi));
tskip=0;
dotSize=30;
figure(1);
PlotPeriodPhasePhaseAsym(EruptionTimes, dotSize, "wtv", 'IDK', false)

figure(3);
clf();

subplot(1,2,1);
sz=40;
scatter(EruptionTimes(:,1), EruptionTimes(:,2), sz, 'filled')

ylim([tskip,Tfinal])


[avgPhase, asymPhase] = eric_phase(EruptionTimes);





indNaN=find(isnan(avgPhase));
indNotNaN=find(~isnan(avgPhase));


cla




subplot(1,3,1);
hold on
% scatter(EruptionTimes(indNotNaN,1),EruptionTimes(indNotNaN,2),dotSize,asymPhase(indNotNaN),'filled','MarkerEdgeColor',[0.8 0.8 0.8])
% scatter(EruptionTimes(indNaN,1),EruptionTimes(indNaN,2),dotSize,[0.6 0.6 0.6])
map = customcolormap_preset('purple-white-green');
map=sqrt(map);
colormap(gca, map);
colorbar

% hold off
caxis(2/3*[-0.3 0.3])
subplot(1,3,2);

hold on
% scatter(EruptionTimes(indNotNaN,1),EruptionTimes(indNotNaN,2),dotSize,asymPhase(indNotNaN),'filled','MarkerEdgeColor',[0.8 0.8 0.8])
% scatter(EruptionTimes(indNaN,1),EruptionTimes(indNaN,2),dotSize,[0.6 0.6 0.6])

map=[127,191,123
    247,247,247;
    175,141,195]/255;

colormap(gca,padmap(map,30))


colorbar
caxis(2/3*[-0.3 0.3])
% hold off


iskip=1+tskip/dt;
tmesh = linspace(0,Tfinal, Nsteps+1);
peaks = get_peak_inds(mod(theta(iskip:end,:), 2*pi),'MinPeakProminence',pi/4);
[ovec,tvec,dphi] = get_phase_asymmetry(peaks, tmesh(tmesh>=tskip),'MinPeakProminence',1);

subplot(1,3,3);

scatter(ovec-(N+1), tvec, dotSize, dphi, "filled");

subtitle('Phase Asymmetry')
xlabel('oscillator #');
ylabel('time')
% xlim([1,N]-mid)
pad = 10;
map=[127,191,123
    247,247,247;
    175,141,195]/255;

colormap(gca,padmap(map,50))

% map = customcolormap_preset('purple-white-green');
% map=sqrt(map);
% colormap(gca, map);
colorbar
clim([-0.2, 0.2])