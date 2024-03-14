function [map,T] = peaktimes_to_cycles(peak_times, tmesh, opts)
arguments
    peak_times (1,:) {mustBeA(peak_times, 'cell')}
    tmesh (1,:) {mustBeNumeric}
    opts.cycle_min (:,:) {mustBeNumeric} = [];
    opts.tmin (1,1) {mustBeNumeric} = -Inf;
    opts.tmax (1,1) {mustBeNumeric} = Inf;
    opts.tburn (1,1) {mustBeNumeric} = 0;
    opts.ind (1,1) {mustBeInteger, mustBePositive}=1;
end

tmin = max(opts.tmin, tmesh(1));
peaks_ref = peak_times{:,opts.ind};

T = mean(diff(peaks_ref(peaks_ref >= opts.tmin & peaks_ref<opts.tmax)));


cycle_min = floor((max(tmin, min([peak_times{:}]))+opts.tburn)/T);
cycle_min = floor(opts.tburn/T);

cycle0 = ceil((peaks_ref(1)-tmin)/T);
cycle_inds = 0:(length(peaks_ref)-1);
map = @(t) cycle_min + interp1(peaks_ref,cycle_inds, t,'linear','extrap');


end