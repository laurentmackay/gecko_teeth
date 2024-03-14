function peaks = get_peak_inds(x, varargin)
N=size(x,2);
peaks=cell(N,1);
for i=1:N
[~, p] = findpeaks(x(:,i),varargin{:});

peaks{i} = p;
end
end