L=1;
N=15;
N_inter = 20;
pad=N_inter/2;

% pdepe
[xmesh, osc_inds] = oscillator_grid(N, N_inter, L, pad=pad);
xosc=xmesh(osc_inds);
dosc = diff(xosc(1:2));
Nx=length(xmesh);
dx=L/(Nx-1);



T=30;

tmesh=linspace(0,150*T,8000);



p=10/T;
k=1;

% space-constant == sqrt(D/p) == k*dosc
D=p*(k*dosc)^2
% D=1e-3


% T=T*(1+2*0.1*abs(xmesh(osc_inds)-L/2)/L);

Imax=N/p;
eps= 0.27*2*pi/(min(T)*Imax)
% eps=0;
IC = zeros(N+Nx,1);
% IC(1:2:N)=pi;
IC(1:N)= randn(N,1)*0.001*2*pi;

% IC=randn(N+Nx,1)*pi*0;
IC(N+1:end)=abs(IC(N+1:end))+randn(Nx,1)*0.01;

J_oo = speye(N,N);
J_of = sparse(1:N,osc_inds', ones(size(osc_inds)), N, Nx);
J_ff = spdiags(ones(Nx,3),-1:1,Nx,Nx);

Jpattern = osc_field_Jpattern(N,N_inter, pad=pad);

options=odeset('AbsTol',1e-4,'RelTol',1e-4, 'MaxStep', min(T)/10, 'Jpattern', Jpattern);


[t,sol] = ode15s(phase_inhibitor_pde(T, D, p, eps, L, N, N_inter, pad=pad), tmesh, IC, options );
%%
figure(1);
iskip=1;
display_oscillations(mod(sol(:,1:N), 2*pi), t, t(iskip), fields={sol(:,N+1:end)});

% figure(3);
% plot(1:Nx, max(sol(:,N+1:end),[],1),'LineWidth',3);
% ylabel('Max of oscillatory solution', 'FontSize',16)
% xlabel('grid point #', 'FontSize',16)
% axis tight
