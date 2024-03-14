function  z = cartesian_to_kuramoto(sol)
N=size(sol,2)/2;
theta=mod(atan2(sol(:,N+1:end),sol(:,1:N)),2*pi);

z=sum(exp(1i*theta),2)/N;
end