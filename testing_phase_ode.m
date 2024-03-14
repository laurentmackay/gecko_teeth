
% xosc=xosc(2:end-1);

N0= 2;
N=2*(N0)+1;
mid=N0+1;



T=2*pi;
tmesh=linspace(0, 4*T,3e3);

eps=2;
m=1;


T=T*(1+2*0.0*abs((0:N-1)/(N-1)-0.5))';

IC=randn(N,1)*0.1*2*pi;
IC(1:2:end-1) = pi/2;
IC(2:2:end) = -pi/2;

options=odeset('AbsTol',1e-9,'RelTol',1e-9, 'MaxStep', min(T)/4);
rhs = phase_inhibitor_nondim_ode(N, 0, eps, m);

tburn=12e4;
tstep=5e4;
ICP=IC;
for i=1:tburn/tstep
[~,sol] = ode15s(rhs, [0, tstep], IC , options);
ICP=IC;
end
[t,sol] = ode15s(rhs, tmesh, sol(end,:) , options);

%%
figure(1);
tskip=0.0*tmesh(end);
iskip = find(tmesh>=tskip,1);
sol_theta = mod(sol(iskip:end,:), 2*pi);
display_oscillations(sol_theta, tmesh(iskip:end), tskip, nrows=3);

return 
%%
figure(3)
plot(t(iskip:end)/(2*pi),sol_theta(iskip:end,1:2), 'LineWidth',2)
legend('\theta_{\rm odd}','\theta_{\rm even}','FontSize',14)
xlabel('t','Interpreter','tex','FontSize',16);
ylabel('\theta','Interpreter','tex','FontSize',16);
set(gca,"YTick",[0,pi, 2*pi], "YTickLabel",["0", "\pi", "2\pi"]);
ylim([0,2*pi])
xlim tight
% exportgraphics(gcf,'anti_phase_ode.pdf')

%%

figure(4)
plot(t/(2*pi),sol_theta, 'LineWidth',2)
legend('\theta_{1}','\theta_{2}','\theta_{3}')
xlabel('t (nondim)','Interpreter','tex');
ylabel('\theta','Interpreter','tex');
axis tight
% exportgraphics(gcf,'spatially_patterned_ode.pdf')
%%
solout = [cos(sol(iskip:end,1:N)), sin(sol(iskip:end,:))];
write_per_sol(sol,tmesh, sol(end,:), rhs , 'auto/ODE/phase/five_phase_anti.dat', NTST=600, options=options);


