N0=40;
N=2*N0+1;
L=N;

N_inter =   20;
pad=N_inter/2;

[xmesh, ~] = oscillator_grid(N, N_inter, L, pad=pad, center='interface');
Nx=length(xmesh);


% D=700;





T=30;
p=0.1;

a = 0.169;
k = 0.6;
b=0.2;
e=1e-2;

f0=0.06;

kappa = .5;
eps = 1;


Jpattern = osc_field_Jpattern(N, N_inter, pad=pad, node=2, center='interface');

IC=abs(randn(2*N+Nx,1)*0.08);
IC(2*N+1:end)=f0+0.06+abs(IC(2*N+1:end))*0.1;
IC0=IC;


% 
sensing_phase = -1;
gamma=0.1;

ode_rhs_and_gamma = FHN_secrete(a, b, e, k, eps, f0=f0, timescale = 0.2*v_gradient(gamma,N), sense_phase=sensing_phase);
rhs = flux_jump_inhibitor_scheme(ode_rhs_and_gamma, kappa, p, N, N_inter,  pad=pad);



% Jac = FHN_grad_compartment_pde_Jac(b, e, 0.1, D, p, L, N, N_inter, pad=pad, grad_dist_func = @(x) 2*abs(x));
% options=odeset('AbsTol',tol,'RelTol',tol, 'MaxStep', min(T)/4, 'Jacobian', Jac);
%%
tol=1e-8;
options=odeset('AbsTol',tol,'RelTol',tol, 'MaxStep', min(T)/4, 'Jpattern', Jpattern, 'Vectorized', 'off');





tic;
ICP=burnin(rhs, IC, tburn=1e-4, tstep=1e4, options=options);
% ICP=IC;
% for i=1:(tburn/tstep)
%     [~, sol] = ode15s(rhs , [0, tstep], ICP, options );
%     ICP=sol(end,:);
% end
toc
%%
tmesh=0:0.5:400*T;



mid=ceil(N/2);
i0=mid+8;
inds = i0:i0+2;
% inds=[];
ablate=false(N,1);
ablate(inds)=true;

% rhs = flux_jump_inhibitor_scheme(ode_rhs_and_gamma,  D, p, L, N, N_inter,  pad=pad, ablate=ablate);



tic 
[t,sol] = ode15s(rhs , tmesh, ICP, options );
toc
IC=sol(end,:);
% ic = cell2mat(arrayfun(bric(a,b,0.01), xmesh, 'UniformOutput', false));
% sol = pbcpdeSolver(brpde_pbc(a,b,Du,Dv), ic, xmesh, tmesh);
%%
figure(3);
tskip = 0.05*tmesh(end);
iskip=find(t>tskip,1);
display_oscillations(sol(:,1:N), t, tskip,checker=true);%, fields={sol(:,(2*N+1):end)});
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

%%
% figure(6);
% plot(sol(:,1:N), sol(:,N+1:2*N)); 
% 
% % hold('on');
% % u=linspace(min(min(sol(:,1:N))),max(max(sol(:,1:N)))); 
% % m=12;
% % step=smooth_heav(); 
% % h=1-step(u,0.1,0.01);
% % vmax=max(max( sol(:,N+1:2*N))); 
% % plot(u,h*vmax/max(h)); 
% % hold('off');
% 
% 
% ylabel('v: inactivation variable ');
% xlabel('u: activation variable ');
% text(0.9,0.08,{'phase 1', '"active"'}, 'Interpreter','tex','HorizontalAlignment','center')
% text(0.27,0.15,'phase 2', 'Interpreter','tex')
% text(-0.15,0.06,{'phase 3', '"inactive"'}, 'Interpreter','tex','HorizontalAlignment','center')
% text(0.6,-0.01,'phase 4', 'Interpreter','tex')
% disp('done')
% return 
%%
usol = sol(iskip:end,1:N);
vsol = sol(iskip:end,N+1:2*N);

umin  = min(usol(:));
umax=max(usol(:));
urange=umax-umin;

vmin  = min(vsol(:));
vmax=max(vsol(:));
vrange=vmax-vmin;

u0 = linspace(umin-0.1*urange,umax+0.1*urange, 100);
v0 = linspace(vmin-0.1*vrange,vmax+0.1*vrange, 100);
[u,v] = meshgrid(u0, v0);

% [x0,n] = FHN_sensing_vectors(sensing_phase, f0=f0);
% u_nullcline = FHN_u_nullcline(a, f0=f0);
% v_nullcline = FHN_v_nullcline(b,k);
% n=n/norm(n);
% step = smooth_heav();
sense =  FHN_sensing_func(sensing_phase, f0=f0);
func = reshape(sense(u(:),v(:)),size(u));
figure(7);
s=pcolor(u,v, func);
set(s,'LineStyle','none')
colorbar
hold on;
plot(u0,u_nullcline(u0),'-c', u0,v_nullcline(u0),'-r');
plot(usol, vsol); 

hold off;

return
%%
peaks = get_sol_peaks(sol(:,1:N));


[ovec,tvec,phi_l, phi_r, phi_ls, phi_rs] = get_phases(peaks, tmesh);
[chi, chis] = phase_asymmetry(phi_l, phi_r, phi_ls, phi_rs);
asymmetry_mean = cellfun(@(p) mean(p), chis);
asymmetry_err = cellfun(@(p) std(p), chis);

[psi, psis] = phase_checker(phi_l, phi_r, phi_ls, phi_rs);
checker_mean = cellfun(@(p) mean(p), psis);
checker_err = cellfun(@(p) std(p), psis);
i_left = ovec<N/2;
i_right = ovec>N/2;
figure(8);
scatter(psi(i_left), chi(i_left),5,'filled', 'Color', 'b' );
hold on
scatter(psi(i_right), chi(i_right),5,'filled', 'Color', 'r' )
hold off
ccl=corrcoef(chi(i_left), psi(i_left));
ccl(2)
ccr=corrcoef(chi(i_right), psi(i_right));
ccr(2)

%%
rectangle('Position',[inds(1)-ceil(N/2) tskip inds(end)-inds(1) tmesh(end)-tskip], 'EdgeColor','r')
%%
get_per_sol(sol(iskip:end,1:N),tmesh, sol(end,:), rhs , NTST=600, options=options);

%%
write_per_sol(sol(iskip:end,1:N),tmesh, sol(end,:), rhs , 'auto/PDE/FHN_secrete/five_6_anti.dat', NTST=600, options=options);
