function theta = bilateral_phase_waves(N,chi, opts)
arguments
N (1,1) {mustBeInteger, mustBePositive}
chi (1,1) {mustBeNumeric}
opts.delta0 (1,1) {mustBeNumeric} = pi;
opts.phi_mid (1,1) {mustBeNumeric} = 0;
end

mid = ceil(N/2);

theta=1:N;

theta_r = phase_waves(mid-1, chi, delta0=opts.delta0);

theta(mid) = opts.phi_mid;
theta(mid+1:end) =  theta_r;
theta(1:mid-1) =  fliplr(theta_r);


end