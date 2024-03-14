L=1;
N=15;
N_inter =   16;
pad=N_inter/2;

[xmesh, osc_inds] = oscillator_grid(N, N_inter, L, pad=pad, center='interface');
xosc=xmesh(osc_inds);
if N>1
dosc = diff(xosc(1:2));
else
    dosc=L/2
end
% D=700;
T=125;

tmesh=linspace(0,12e3,12e3);
% (1+cos(theta))/2*(N)/p
p=0.01;
% p=0;
kappa=3;
D=p*(kappa*dosc)^2
% D=0.0008
% D=0.005
% D=4e-6;
% D=1e-4;

Nx=length(xmesh);

dx=L/(Nx-1);

x=linspace(0,L,Nx)';


% T=T*(1+2*0.1*abs(x(osc_inds)-L/2)/L);
e=1e-2;

% e=e*(1-2*0.1*abs(x(osc_inds)-L/2)/L);

Jpattern = osc_field_Jpattern(N, N_inter, pad=pad, node=2, center='interface');

IC=randn(2*N+Nx,1)*0.8;
IC(2*N+1:end)=abs(IC(2*N+1:end))*0.001;
for i=1:N
Jpattern(i,2*N+1:end)=0;
Jpattern(2*N+1:end,i)=0;
end
options=odeset('AbsTol',1e-8,'RelTol',1e-8, 'MaxStep', min(T)/4, 'JPattern',Jpattern);
% options=odeset('AbsTol',1e-8,'RelTol',1e-8, 'MaxStep', min(T)/4);


b=0.2;
% D=35e-9;
tic;
rhs=SCH_grad_compartment_pde(b,  -0.1, D, p, L, N, N_inter, pad=pad, grad_dist_func = @(x) 2*abs(x));

[~,sol] = ode15s(rhs , [0, 1e4], IC, options );

    %%
tmesh=linspace(0,4e3,16e4);

[t,sol] = ode15s(rhs , tmesh, sol(end,:), options );
toc
% ic = cell2mat(arrayfun(bric(a,b,0.01), xmesh, 'UniformOutput', false));
% sol = pbcpdeSolver(brpde_pbc(a,b,Du,Dv), ic, xmesh, tmesh);
%%
figure(4);
iskip=4e4;
display_oscillations(sol(:,1:N), t, t(iskip), fields={sol(:,(2*N+1):end)});
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
return 

%%
write_per_sol(sol(iskip:end,:),tmesh, sol(end,:), rhs , 'three_FHN_anti.dat', NTST=600, options=options);
