L=1;
N=48;
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
dx=xmesh(2)-xmesh(1)
tmesh=linspace(0,12e3,12e3);
% (1+cos(theta))/2*(N)/p
p=0.01;
% p=0;
kappa=1.5;
D=p*(kappa*dosc)^2
% D=0.0008
% D=4e-6;
% D=1e-4;

Nx=length(xmesh);

dx=L/(Nx-1);

x=linspace(0,L,Nx)';


% T=T*(1+2*0.1*abs(x(osc_inds)-L/2)/L);
e=1e-2;

% e=e*(1-2*0.1*abs(x(osc_inds)-L/2)/L);

Jpattern = osc_field_Jpattern(N, N_inter, pad=pad, node=2, center='interface');

IC=randn(2*N+Nx,1)*0.08;
IC(2*N+1:end)=abs(IC(2*N+1:end))*0.001;
% for i=1:N 
% Jpattern(i,2*N+1:end)=0;
% Jpattern(2*N+1:end,i)=0;
% end



% options=odeset('AbsTol',1e-8,'RelTol',1e-8, 'MaxStep', min(T)/4);


b=.1;
% D=35e-9;

rhs=FHN_grad_compartment_pde(b, e, 0.1, D, p, L, N, N_inter, f0=0.08, pad=pad, grad_dist_func = @(x) x);


% Jac = FHN_grad_compartment_pde_Jac(b, e, 0.1, D, p, L, N, N_inter, pad=pad, grad_dist_func = @(x) 2*abs(x));
% options=odeset('AbsTol',tol,'RelTol',tol, 'MaxStep', min(T)/4, 'Jacobian', Jac);

tol=1e-5;
options=odeset('AbsTol',tol,'RelTol',tol, 'MaxStep', min(T)/4, 'Jpattern', Jpattern, 'Vectorized', 'off');





tic;
ICP=burnin(rhs, IC, tburn=4e-4, tstep=1e-4, options=options);
% ICP=IC;
% for i=1:(tburn/tstep)
%     [~, sol] = ode15s(rhs , [0, tstep], ICP, options );
%     ICP=sol(end,:);
% end
toc
    %%
tmesh=0:1:30e3;


i0=ceil(N/2) + 6;
inds = i0:(i0+5);
% ICP(inds)=0;
% ICP(N+inds)=0;


[t,sol] = ode15s(rhs , tmesh, ICP, options );

% ic = cell2mat(arrayfun(bric(a,b,0.01), xmesh, 'UniformOutput', false));
% sol = pbcpdeSolver(brpde_pbc(a,b,Du,Dv), ic, xmesh, tmesh);
%%
figure(4);
tskip = 0.01*tmesh(end);
iskip=find(t>tskip,1);
display_oscillations(sol(:,1:N), t, tskip , checker=true);%, fields={sol(:,(2*N+1):end)});
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
rectangle('Position',[inds(1)-ceil(N/2) tskip inds(end)-inds(1) tmesh(end)-tskip], 'EdgeColor','r')
%%
write_per_sol(sol(iskip:end,:),tmesh, sol(end,:), rhs , 'three_FHN_anti.dat', NTST=600, options=options);
