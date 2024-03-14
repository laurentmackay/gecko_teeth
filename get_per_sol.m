function [sol,T] = get_per_sol(samp, t, IC, rhsfun, opts)
arguments
    samp  (:,:) {mustBeNumeric} 
    t  (:,1) {mustBeNumeric} 
    IC  (:,1) {mustBeNumeric} 
    rhsfun (1,1) function_handle
    opts.options (1,1) struct = odeset()
    opts.NTST (1,1) {mustBeNumeric}  = 150
    opts.periodicity  (1,1) {mustBeInteger} = 1
    opts.MinPeakProminence (1,1) {mustBeNumeric} = 0.1
end
figure(2)
Nsamp = size(samp,2);
per = zeros(Nsamp,1);
for i=1:Nsamp
    vmax =  max(samp(:,i));
    vmin =  min(samp(:,i));
    range = vmax-vmin;
    findpeaks(samp(:,i),'MinPeakProminence',opts.MinPeakProminence*range)
[~,p] = findpeaks(samp(:,i),'MinPeakProminence',opts.MinPeakProminence*range);
per(i) = mean(diff(t(p-1)));
drawnow
end
per

T=mode(per)*opts.periodicity

T=mean(per)*opts.periodicity

tspan = linspace(0,T, opts.NTST);
% options = odeset('MaxStep',T);
[T, sol] = ode15s( rhsfun ,tspan, IC, opts.options);

rel_delta = (sol(1,:) - sol(end,:))./ max(abs(sol),[],1);
disp(['max rel. delta = ' num2str(max(abs(rel_delta)))])
sol(end,:)=sol(1,:);
end