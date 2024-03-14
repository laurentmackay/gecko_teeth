function theta = phase_asym_IC3(N,chi, opts)
arguments
N (1,1) {mustBeInteger, mustBePositive}
chi (1,1) {mustBeNumeric}
opts.delta0 (1,1) {mustBeNumeric} = pi;
end

theta=1:N;

theta(1:2:N) = cumsum(repmat(chi*2*pi,[ceil(N/2),1]))+opts.delta0;
theta(2:2:N) = cumsum(repmat(chi*2*pi,[floor(N/2),1]));
theta = mod(theta, 2*pi);
end