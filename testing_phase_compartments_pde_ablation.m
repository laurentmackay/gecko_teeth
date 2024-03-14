L=1;
N=31;
N_inter =   12;
pad=6*N_inter;

[xmesh, osc_inds] = oscillator_grid(N, N_inter, L, pad=pad, center='interface');

Nx=length(xmesh);
xosc=(xmesh(osc_inds)+xmesh(osc_inds+1))/2;

if N>1
dosc = diff(xosc(1:2));
else
    dosc=L/2
end





T0 = 30;
p=.1/T0;

kappa=2;
D=p*(kappa*dosc)^2





Iave = 0.5*N/p;
eps=.25*2*pi/(Iave*T0)
% D=35e-9;


IC=randn(2*N+Nx,1)*0.0001;
odd = [1:2:N, N+1:2:2*N];
even = [2:2:N-1, N+2:2:2*N-1];
IC(odd)=sqrt(2)/2;
IC(even)=-sqrt(2)/2;
IC(1:2*N) = IC(1:2*N) + 0.001*randn(2*N,1);
theta0 = atan2(IC(N+1:2*N),IC(1:N));
% IC(2*N+1:end)
% IC(2*N+1:3*N) = secrete(theta0);
IC(2*N+1:end)=Iave+IC(2*N+1:end);
% IC(2*N+1:end)=abs(IC(2*N+1:end))*0.001;

Jpattern = osc_field_Jpattern(N, N_inter, pad=pad, node=2, center='interface');
% for i=1:N
%     Jpattern(i,2*N+1:end)=0;
%     Jpattern(2*N+1:end,i)=0;
% end

tol=1e-5;
options=odeset('AbsTol',tol,'RelTol',tol, 'MaxStep', min(T0)/4, 'JPattern',Jpattern);
% options=odeset('AbsTol',1e-8,'RelTol',1e-8, 'MaxStep', min(T)/4);

gamma = 0.001;


rhs=cartesian_phase_compartment_pde(eps, T0, gamma, D, p, L, N, N_inter, pad=pad, grad_dist_func = @(x) 2*abs(x));


%%
tic;
ICP=burnin(rhs, IC, tburn=30e4, tstep=1e4, options=options);
    %%
% 
% ICP = ICP_lt;

% ICP(inds)=1;
% ICP(N+inds)=0;
mid=ceil(N/2);
i0=mid+15;
inds = i0:i0+5;
inds=[];
ablate=false(N,1);
ablate(inds)=true;
 % load('small_grad.mat')
rhs=cartesian_phase_compartment_pde(eps, T0, gamma, D, p, L, N, N_inter, pad=pad, grad_dist_func = @(x) 2*abs(x), ablate=ablate);


tmesh=0:1:600*T0;

[t,sol] = ode15s(rhs , tmesh, ICP, options );
% ICP=sol(end,:);
toc
% ic = cell2mat(arrayfun(bric(a,b,0.01), xmesh, 'UniformOutput', false));
% sol = pbcpdeSolver(brpde_pbc(a,b,Du,Dv), ic, xmesh, tmesh);
%%
figure(5);
tskip=0.00*t(end);
iskip=find(t>=tskip,1);
sol_theta=sol;
sol_theta(:,1:N) = 0.5+atan2(sol(:,N+1:2*N),sol(:,1:N))/(2*pi);
display_oscillations(sol_theta(:,1:N), t, tskip, fields = {sol(:,2*N+1:end)});
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
if ~isempty(inds)
    rectangle('Position',[inds(1)-ceil(N/2)-0.5 tskip inds(end)-inds(1)+1 tmesh(end)-tskip], 'EdgeColor','r');
end
disp('done')
return 
%%

%%
get_per_sol(sol(iskip:end,1:2*N), tmesh, sol(end,:), rhs, NTST=600, options=options);
%%
write_per_sol(sol(iskip:end,1:2*N), tmesh, sol(end,:), rhs, ['auto/PDE/phase/seven_anti_gradient.dat'], NTST=600, options=options);
 