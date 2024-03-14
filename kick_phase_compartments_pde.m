
N0=5;
N=2*N0+1;
L=N;
N_inter =   20;
pad=N_inter/2;

[xmesh, osc_inds] = oscillator_grid(N, N_inter, L, pad=pad, center='interface');

Nx=length(xmesh);
xosc=xmesh(osc_inds);







T0 = 30;
p=1/T0;

kappa=5;

% kappa=1.5;





f=.95;
% f=0.7;
% eps=f*2*pi/(T0*Iave)
% D=35e-9;


IC=randn(2*N+Nx,1)*0.0001;
odd = [1:2:N; N+1:2:2*N];
even = [2:2:N-1; N+2:2:2*N-1];
t0 = pi/4 + 0.00*randn(N,1);
% mid=ceil(N/2)
mid=ceil(N/2);
adv=0.25;
t0(N0+2:end) =  phase_asym_IC(N0, 1-adv, 0.2);
t0(1:N0) =  phase_asym_IC(N0, adv, -0.2);

IC(odd(1,:))=cos(t0(odd(1,:)));
IC(odd(2,:))=sin(t0(odd(1,:)));
IC(even(1,:))=cos(t0(even(1,:)));
IC(even(2,:))=sin(t0(even(1,:)));
% IC(even)=sin(pi+t0);
% IC(1:2*N) = IC(1:2*N) 
theta0 = atan2(IC(N+1:2*N),IC(1:N));
% IC(2*N+1:end)
% IC(2*N+1:3*N) = secrete(thet2a0);
IC(2*N+1:end)=1+IC(2*N+1:end);
% IC(2*N+1:end)=abs(IC(2*N+1:end))*0.001;

Jpattern = osc_field_Jpattern(N, N_inter, pad=pad, node=2, center='interface');
% for i=1:N
%     Jpattern(i,2*N+1:end)=0;
%     Jpattern(2*N+1:end,i)=0;
% end

tol=1e-10;
options=odeset('AbsTol',tol,'RelTol',tol, 'MaxStep', min(T0)/4, 'JPattern',Jpattern);
% options=odeset('AbsTol',1e-8,'RelTol',1e-8, 'MaxStep', min(T)/4);

% 
% gamma = -0.2;
% rhs=cartesian_phase_compartment_pde_pgrad(f, T0, gamma , kappa, p, N, N_inter, pad=pad, grad_dist_func = @(x) 2*abs(x));
gamma=0.0;
rhs=cartesian_phase_compartment_pde_nondim(f, T0, gamma , kappa, p, N, N_inter, pad=pad, grad_dist_func = @(x) 2*abs(x));


%%
tic;
ICP=burnin(rhs, IC, tburn=4e4, tstep=4e3, options=options);
toc
    %%
theta0 = atan2(ICP(N+1:2*N),ICP(1:N));

kick=0.25*2*pi;



ICP(1)=cos(theta0(1)+kick);
ICP(N+1)=sin(theta0(1)+kick);

tmesh=0:.1:145*T0;
tic;
[t,sol] = ode15s(rhs , tmesh, ICP, options );
ICP=sol(end,:);
IC=ICP;
toc
% ic = cell2mat(arrayfun(bric(a,b,0.01), xmesh, 'UniformOutput', false));
% sol = pbcpdeSolver(brpde_pbc(a,b,Du,Dv), ic, xmesh, tmesh);
%%
figure(1);
tskip=0*t(end);
iskip=find(t>=tskip,1);
% sol_theta=sol;
sol_theta = atan2(sol(:,N+1:2*N),sol(:,1:N));
display_oscillations(sol_theta(:,1:N), t, tskip, checker=false, fields = {sol(:,2*N+1:end)});
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
R=5;
r=2;
imax=8e3;
theta = sol_theta(1:imax,1);
phi = sol_theta(1:imax,N0);
figure(15);plot3((R+r*cos(theta)).*cos(phi), (R+r*cos(theta)).*sin(phi), r*sin(theta))
%%
 figure(15);plot3(sol(:,2*N+osc_inds(N0+5)), sol(:,2*N+osc_inds(N0)), sol(:,2*N+osc_inds(N0+5)-1))
%%
disp(u_str(sol(end,:)))

%%
get_per_sol(sol(iskip:end,1:2*N), tmesh, sol(end,:), rhs, NTST=600, options=options);
%%
write_per_sol(sol(iskip:end,1:2*N), tmesh, sol(end,:), rhs, ['auto/PDE/phase/five_kappa=0.5_f=0.4_gamma=0.002_sync.dat'], NTST=600, options=options);
 