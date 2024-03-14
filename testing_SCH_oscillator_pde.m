L=1;
N=15;
N_inter = 10;
pad=N_inter/2;

[xmesh, osc_inds] = oscillator_grid(N, N_inter, L, pad=pad);
xosc=xmesh(osc_inds);
dosc = diff(xosc(1:2));

% D=700;
T=300;
 dx=xmesh(2)-xmesh(1)
tmesh=linspace(0,4e3,1e3);
% (1+cos(theta))/2*(N)/p
p=0.1;
% p=0;
k=2;
D=p*(k*dosc)^2
% D=1
% D=4e-6;


Nx=length(xmesh);

dx=L/(Nx-1);

x=linspace(0,L,Nx)';


% T=T*(1+2*0.1*abs(x(osc_inds)-L/2)/L);
e=1e-2;
% e=e*(1-2*0.1*abs(x(osc_inds)-L/2)/L);

Jpattern = osc_field_Jpattern(N, N_inter, pad=pad);

IC=randn(N+Nx,1)*0.08;
IC(N+1:end)=abs(IC(N+1:end))*0.001;
options=odeset('AbsTol',1e-5,'RelTol',1e-5, 'MaxStep', min(T)/4,'JPattern',Jpattern);

b=.2;
% D=35e-9;

[t,sol] = ode15s(SCH_oscillator_pde(b, D, p, L, N, N_inter, pad=pad), tmesh, IC, options );
% ic = cell2mat(arrayfun(bric(a,b,0.01), xmesh, 'UniformOutput', false));
% sol = pbcpdeSolver(brpde_pbc(a,b,Du,Dv), ic, xmesh, tmesh);
%%
figure(1);
iskip=1;
display_oscillations(sol(:,1:N), t, t(iskip), fields={sol(:,(N+1):end)});
% subplot(1,2,1);
% s=imagesc(1:N, flipud(t(iskip:end)),mod(sol(iskip:end,1:N), 2*pi));
% set(gca,'YTickLabel',flipud(get(gca, 'YTickLabel')))
% % s.EdgeAlpha=0;
% colorbar
% subplot(1,2,2);
% s=imagesc(xmesh, t(iskip:end),flipud(eps*sol(iskip:end,N+1:end))/(2*pi/T));
% set(gca,'YTickLabel',flipud(get(gca, 'YTickLabel')))
% % s=pcolor(xmesh, tmesh(iskip:end), sol(iskip:end,:,2));
% % s.EdgeAlpha=0;

disp('done')