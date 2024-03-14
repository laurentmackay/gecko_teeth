
N0=22;
N=2*N0+1;
L=N;
N_inter =  20;
pad=floor(N_inter/2);

[xmesh, osc_inds] = oscillator_grid(N, N_inter, L, pad=pad, center='interface');

Nx=length(xmesh);
xosc=xmesh(osc_inds);






T0 = 30;

omega=1;
kappa=1;
f=0.1;
gamma=0;

tburn = 2e4*T0;


p = 2*pi/(omega*T0);
msense=1;
msecrete=1;
% f=0.775;
% eps=f*2*pi/(T0*Iave)
% D=35e-9;


IC=randn(2*N+Nx,1)*0.01;
odd = [1:2:N; N+1:2:2*N];
even = [2:2:N-1; N+2:2:2*N-1];
t0 = pi/4 + 0.00*randn(N,1);
% mid=ceil(N/2)
mid=ceil(N/2);
adv=0.15;
chi = 0.2;
flip=1;
mid

theta_r = phase_asym_IC3(N0, 0.1);
t0(N0+1) = 0;
t0(N0+2:end) =  theta_r;
t0(1:N0) =  fliplr(theta_r);
% 
% t0(N0+2:N0+2+floor(N0/2-1)) =  fliplr(theta_r(1:1+floor(N0/2-1)));
% t0(N0-floor(N0/2-1):N0) =  theta_r(1:1+floor(N0/2-1));


% t0(mid+1:end) =  phase_asym_IC(N0, -flip*adv, flip*chi);
% t0(1:N0) =  phase_asym_IC(N0, flip*adv, -flip*chi);
% t0(mid+1:end) =  phase_asym_IC(N0, -flip*adv, flip*chi);
% t0(:) =  phase_asym_IC(N, flip*adv, -flip*chi);
x0=0;

y0=0;

x0=cos(t0(mid-1))+cos(t0(mid+1));
y0=sin(t0(mid-1))+sin(t0(mid+1));
% t0(mid) = atan2(y0,x0)+pi;
t0 = t0 + 0*randn(size(t0));

% thetal= mod(atan2(ICl(N+1:2*N),ICl(:,1:N)),2*pi);
% thetar= mod(atan2(ICr(N+1:2*N),ICr(:,1:N)),2*pi);
% 
% t0(1:mid-1)= thetal(1:mid-1);
% t0(mid+1:end)=thetar(mid+1:end);
% t0(mid)=mod(pi+(thetar(mid+1)+thetal(mid-1))/2,2*pi);
% 


% t0 = pi/2 + 0.1*randn(N,1);

% t0 =  phase_asym_IC(N, adv, 0.2)+ 0.01*randn(N,1);
% % 
% t0(odd(1,:)) = t0(odd(1,:)) + pi ;
% % 



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


% IC(2*N+1:osc_inds(mid))=ICl(2*N+1:osc_inds(mid));
% % IC(2*N+osc_inds(mid):end)=ICr(2*N+osc_inds(mid):end);
% IC=(ICl+ICr)/2;
% tburn=0.0001;
% IC=IC_acc;


Jpattern = osc_field_Jpattern(N, N_inter, pad=pad, node=2, center='interface');


tol=1e-10;
options=odeset('AbsTol',tol,'RelTol',tol, 'MaxStep', min(T0)/4, 'JPattern',Jpattern);
% options=odeset('AbsTol',1e-8,'RelTol',1e-8, 'MaxStep', min(T)/4);

% 
% gamma = -0.2;
% rhs=cartesian_phase_compartment_pde_pgrad(f, T0, gamma , kappa, p, N, N_inter, pad=pad, grad_dist_func = @(x) 2*abs(x));


% rhs0=cartesian_phase_compartment_pde_nondim(f, 1e15, 0 , kappa, p, N, N_inter, pad=pad, grad_dist_func = @(x) 2*abs(x), msense=msense, msecrete=msecrete);


% IC=burnin(rhs0, IC, tburn=6e3, tstep=2e3, options=options);
% 
% IC0=IC;
mper=1

rhs=cartesian_phase_compartment_pde_nondim(f, T0, gamma , kappa, p, N, N_inter, pad=pad, grad_dist_func = @(x) (abs(x)/0.5).^mper, msense=msense, msecrete=msecrete);

% IC=sol0(end,:);

%%
tic;
ICP=burnin(rhs, IC, tburn=tburn, tstep=1e3, options=options);
toc
    %%
% tmesh=linspace(0,2e3,2e5);

mid=ceil(N/2);
i0=mid+15;
inds = i0:i0+6;
% inds=[];
% ablate=false(N,1);
% ablate(inds)=true;
% 
% 
% 
% rhs=cartesian_phase_compartment_pde_nondim(f, T0, gamma, kappa, p, N , N_inter, pad=pad, grad_dist_func = @(x) (abs(x)/0.5).^mper,  msense=msense, msecrete=msecrete);
% 
% 


dt=0.01;
tmesh=0:dt:30*T0;
tic;
[t,sol] = ode15s(rhs , tmesh, ICP, options );
ICP=sol(end,:);
IC=ICP;
toc
% ic = cell2mat(arrayfun(bric(a,b,0.01), xmesh, 'UniformOutput', false));
% sol = pbcpdeSolver(brpde_pbc(a,b,Du,Dv), ic, xmesh, tmesh);
%%

tskip=0.0*t(end);
iskip=find(t>=tskip,1);
tlast=1*t(end);
ilast=find(t>=tlast);
% sol_theta=sol;
sol_theta = mod(atan2(sol(:,(N+1):2*N),sol(:,1:N)), 2*pi);
sol_radius = sqrt(sol(:,(N+1):2*N).^2+sol(:,1:N).^2);
% figure(3)
% display_oscillations(sol_theta(1:ilast,1:N), t(1:ilast), tskip,checker=true, fields = {sol(:,2*N+1:end)});
f0=0.00;
w=.2;
figure(1);
oscillations_summary(sol_theta(1:ilast,1:N), t(1:ilast), tmin=f0*t(end), tmax=(f0+w)*t(end), tburn=-80,...
    xlines=[ 0], border=true, dotSize=25, phi_scatter=[], checker=false, edmund_lines=true, period=false, ...
    cycle_min=0, cycle_max=3);
% figure(2);
% oscillations_summary(sol_theta(1:ilast,1:N), t(1:ilast), tmin=f0*t(end), tmax=(f0+w)*t(end), xlines=[-20:20], border=false, phi_scatter=[], checker=true);


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
save(['croc_forwards_f=' num2str(f) '_omega=' num2str(omega) '_gamma=' num2str(gamma) '_kappa=' num2str(kappa) '_N=' num2str(N) '.mat'],'IC_acc')

%%
fn = ['auto/PDE/phase/images/pos_f=' num2str(f) '_omega=' num2str(omega) '_gamma=' num2str(gamma) '_kappa=' num2str(kappa) '_N=' num2str(N) '.pdf'];
exportgraphics(gcf, fn)
%%
figure(11);
peak_inds = get_peak_inds(sol_theta(iskip:end,:));
locations = -N0:N0;
EruptionTimes = arrayfun(@(i) [repmat(locations(i), [ length(peak_inds{i}),1]), peak_inds{i}], 1:length(peak_inds), 'UniformOutput', false);
EruptionTimes =  vertcat(EruptionTimes{:});
EruptionTimes = sortrows(EruptionTimes,2);
% PlotPeriodPhasePhaseAsym(EruptionTimes, dotSize, "wtv", 'Some Sort of Tile', false)
PlotPeriodPhaseAsym(EruptionTimes, 55, ['pde_model_kappa=' num2str(kappa) '_omega=' num2str(omega) '_epsilon=' num2str(f)], '', false, true)

%%
R=5;
r=2;
imax=8e4;
theta = sol_theta(1:imax,1);
phi = sol_theta(1:imax,N0);
figure(15);plot3((R+r*cos(theta)).*cos(phi), (R+r*cos(theta)).*sin(phi), r*sin(theta))
%%
[~,eLag,eDim] = phaseSpaceReconstruction(sol_theta)
lmax = lyapunovExponent(sol_theta,dt,'ExpansionRange', [1,3000])
%%
figure(34);plot(sol(iskip:end,i), sol(iskip:end,j))
xlabel('I(x_1, t)')
ylabel('I(0, t)','Rotation',0.5)
xbnd = xlim();
ybnd = ylim();
% text(xbnd(1) +(diff(xbnd))*.7, ybnd(1)+diff(ybnd)*.9,['\lambda_{max} = ' num2str(lmax)])
%%

i=2*N+osc_inds(1);
j=2*N+1
k=2*N+osc_inds(N);
figure(15);plot3(sol(iskip:end,i), sol(iskip:end,j), sol(iskip:end,k));
xlabel(['I(x_{' num2str(1) '},t)'])
ylabel(['$$I(\frac{x_1 +x_2}{2},t)$$'],'Interpreter','latex')
zlabel(['$$I(x_{' num2str(N) '},t)$$'],'Interpreter','latex')
xbnd = xlim();
ybnd = ylim();
zbnd = zlim();
text(xbnd(1) +(diff(xbnd))*.25, ybnd(1)+diff(ybnd)*.95, zbnd(1)+diff(zbnd)*.9,['\lambda_{max} = ' num2str(lmax)])

%%
I_osc=sol(:,2*N+osc_inds);
figure(55);plot(mod(sol_theta, 2*pi), I_osc,'.'); hold on; theta=linspace(0,2*pi,100); plot(theta,1+(1+cos(theta))/2); hold off;
%%
disp(u_str(sol(end,:)))

%%
get_per_sol(sol(iskip:end,1:2*N), tmesh, sol(end,:), rhs, NTST=1000, options=options, periodicity=1, MinPeakProminence=0.95);
%%
fn = ['auto/PDE/phase/' num2str(N) '/ftb_kappa=' num2str(kappa) '_f=' num2str(f) '_gamma=' num2str(gamma) '_omega=' num2str(omega) '.dat']
write_per_sol(sol(iskip:end,1:2*N), tmesh, sol(end,:), rhs, fn, NTST=1000, options=options, MinPeakProminence=0.95);
 
%%
exportgraphics(gcf,['double_waves_N=' num2str(N) '_gamma=' num2str(gamma) '_eps=' num2str(f) '_kappa=' num2str(kappa) '_omega=' num2str(omega) '.pdf'])

%%

IC0=IC;
t0p=sol_theta(iskip,:);
t0p(:)=(pi/2+0.1*(t0p-pi/2));
IC0(odd(1,:))=cos(t0p(odd(1,:)));
IC0(odd(2,:))=sin(t0p(odd(1,:)));
IC0(even(1,:))=cos(t0p(even(1,:)));
IC0(even(2,:))=sin(t0p(even(1,:)));
rhs0=cartesian_phase_compartment_pde_nondim(f, 1e15, 0 , kappa, p, N, N_inter, pad=pad, grad_dist_func = @(x) 2*abs(x), msense=msense, msecrete=msecrete);
IClag=burnin(rhs0, IC0, tburn=6e3, tstep=2e3, options=options);

figure(55); plot([sol(iskip,2*N+1:end); IClag(2*N+1:end)]')

