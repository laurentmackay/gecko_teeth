function  display_oscillations(sol, tmesh, tskip, extra)
arguments
    sol (:,:) {mustBeNumeric}
    tmesh (1,:) {mustBeNumeric}
    tskip (1,1) {mustBeNumeric}
    extra.fields (1,:) cell = []
    extra.nrows  (1,1) {mustBeInteger} = 3
end

upper_style = {'XTick',[], 'FontSize',14};
lower_style = { 'FontSize',18};

iskip=find(tmesh>=tskip,1);
sol=sol(iskip:end,:);
tmesh=tmesh(iskip:end);

N=size(sol,2);
if N>1
    mid=ceil(N/2)+(1-mod(N,2))/2;
    oscs = (1:N)-mid;
else
    mid=1;
    oscs=[0];
end



smax=max(sol);
smin=min(sol);
range=smax-smin;

peaks = cell(N,1);
for i=1:N
[~,p]=findpeaks(sol(:,i),'MinPeakProminence',range(i)/4 );
peaks{i} = p;
end 

peak_times = cellfun(@(p) tmesh(p), peaks,'UniformOutput', false);
mean_period = cellfun(@(t) mean(diff(t)), peak_times);
period_err = cellfun(@(t) std(diff(t)), peak_times);




[ovec, tvec, dphi, dphis] = get_phase_asymmetry(peaks, tmesh);
asymmetry_mean = cellfun(@(p) mean(p), dphis);
asymmetry_err = cellfun(@(p) std(p), dphis);
plot(dphis{2},'-o')
hold on
cellfun(@(d)plot(d),dphis(3:end-1))
hold off

amax=cellfun(@(d)max(d),dphis(2:end-1))
amin=cellfun(@(d)min(d),dphis(2:end-1))
arange=amax-amin;

apeaks = cell(N-2,1);
for i=2:N-1
[~,p]=findpeaks(-dphis{i},'MinPeakProminence',arange(i-1)/4 );
apeaks{i-1} = p;
end 

ints = cellfun(@(p) diff(p), apeaks,  'UniformOutput', false)

end