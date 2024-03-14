function h = phase_inhibitor_pde(T, eps, m)


omega_dot_0 = (2*pi)./T;
function dydt = inner(t,y)

yr=y(3:end);
yl=y(1:end-2);

dydt = omega_dot_0 - eps*[secrete(y(2)); (secrete(yr)+secrete(yl)); secrete(y(end-1))].*sense(y);

end


h=@inner;
end