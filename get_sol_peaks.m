function peaks = get_sol_peaks(sol, x, opts)
arguments
    sol  (:,:) {mustBeNumeric} 
    x (1,:) {mustBeNumeric}
    opts.MinPeakDistance (1,1)  {mustBeNumeric}  = 0
    opts.MinPeakProminence (1,1)  {mustBeNumeric}  = 0.1
    opts.MinPeakHeight (1,1)  {mustBeNumeric}  = 1e-5;
    
end

N=size(sol,2);


smax=max(sol);
smin=min(sol);
range=smax-smin;

peaks = cell(N,1);
for i=1:N
    if any(sol(:,i)>opts.MinPeakHeight)
        [~,p]=findpeaks(sol(:,i), x,'MinPeakProminence',range(i)*opts.MinPeakProminence, 'MinPeakDistance', opts.MinPeakDistance , 'MinPeakHeight',opts.MinPeakHeight);
        peaks{i} = p;
    else
        peaks{i}=[];
    end
end 
end