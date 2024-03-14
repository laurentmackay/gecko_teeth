function h = FHN_oscillator_pde(b,  D, p, L, N, N_inter, opts)
arguments
    b  (1,1) {mustBeNumeric} 
    D  (1,1) {mustBeNumeric} 
    p (1,1) {mustBeNumeric} 
    L  (1,1) {mustBeNumeric} 
    N  (1,1) {mustBeInteger} 
    N_inter  (1,1) {mustBeInteger}
    opts.BC (1,:) string = 'Neumann'
    opts.pad (1,1) {mustBeInteger} = N_inter/2
end
a = 0.1;

% e = 1e-2;



% Nx=N+(N+1)*N_inter;


[xmesh, osc_inds] = oscillator_grid(N, N_inter, L, pad=opts.pad);
Nx=length(xmesh);

dx=L/(Nx-1);

osc_inds_abs = osc_inds+N;


function dydt = inner(t,y)

u=y(1:N);
v=y(N+1:end);

flux_r=[-D*diff(v)/dx; 0];
% flux_r(1) = flux_r(1)*2;

flux_l=[0; flux_r(1:end-1)];

v_osc = v(osc_inds);
rx=u.*u.*v_osc;

f_osc = a-u+rx;
g_osc =  b-rx;

% g_osc(:)=0;
% flux_r(osc_inds)=flux_l(osc_inds) - g_osc*dx; %this is a result of integrating around the singular source

div_flux = (flux_r-flux_l)/dx;
div_flux(1) = div_flux(1)*2;
div_flux(end) = div_flux(end)*2;


dydt = [u.*(1-u).*(u-a)-v_osc;   (-p*v - div_flux)];

dydt(osc_inds_abs) =  dydt(osc_inds_abs) + g_osc/dx;

% dydt(osc_inds_abs) = -p*v_osc +  g_osc;

end


h=@inner;
end