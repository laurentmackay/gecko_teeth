function h = phase_inhibitor_pde(N, gamma, eps, m)
i=1:N;
locations = 2*((i-1)/(N-1) - 0.5);
T=1+gamma*abs(locations');
omega_dot_0 = 1./T;
function dydt = inner(t,y)

yr=y(3:end);
yl=y(1:end-2);

dydt = omega_dot_0 - eps*[secrete(y(2)); (secrete(yr)+secrete(yl))/2; secrete(y(end-1))].*sense(y);

end


h=@inner;
end