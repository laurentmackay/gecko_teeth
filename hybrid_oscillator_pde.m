function h = FHN_oscillator_pde(f,g  D, p, L, N, N_inter)



Nx=N+(N-1)*N_inter;




dx=L/(Nx-1);
h_cell = repmat(dx,Nx,1);
h_cell(1) = dx/2;
h_cell(end) = dx/2;
osc_inds = (1:(N_inter+1):Nx)';
osc_inds_abs = osc_inds+N;
% x=linspace(0,L,Nx)';
% 
h_osc=h_cell(osc_inds);

function dydt = inner(t,y)

u=y(1:N);
v=y(N+1:end);

flux_r=[diff(v)/dx; 0];
% flux_r(1) = flux_r(1)*2;

flux_l=[0; flux_r(1:end-1)];
v_osc = v(osc_inds);
flux_osc =  e.*(k*u-v_osc-b);

% flux_r(osc_inds)=flux_l(osc_inds)+h_osc.*flux_osc/D;

d2Idx2 = (flux_r-flux_l)/dx;
d2Idx2(1) = d2Idx2(1)*2;
d2Idx2(end) = d2Idx2(end)*2;


dydt = [f(u,v_osc);   -p*v + D*d2Idx2];
dydt(osc_inds_abs) =  g(u,v_osc);
% dydt(osc_inds_abs) = dydt(osc_inds_abs) + 0.5*flux_osc./h_osc;
end


h=@inner;
end