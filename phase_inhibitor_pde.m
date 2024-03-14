function h = phase_inhibitor_pde(T, D, p, eps, L, N, N_inter, opts)
arguments
    T  (:,1) {mustBeNumeric} 
    D  (1,1) {mustBeNumeric} 
    p  (1,1) {mustBeNumeric} 
    eps  (1,1) {mustBeNumeric} 
    L  (1,1) {mustBeNumeric} 
    N  (1,1) {mustBeInteger} 
    N_inter  (1,1) {mustBeInteger}
    opts.BC (1,:) string = 'Neumann'
    opts.pad (1,1) {mustBeInteger} = N_inter/2
end


switch opts.BC
    case 'Neumann'
        pl=0;
        pr=0;
    case 'Robin'
        pl=-p;
        pr=p;
    otherwise
        disp("Unrecognized boundary condition, allowed options: 'Neumann' and 'Robin'")
end
[xmesh, osc_inds] = oscillator_grid(N, N_inter, L, pad=opts.pad);
Nx=length(xmesh);

dx=L/(Nx-1);

osc_inds_abs = osc_inds+N;

h_cell = repmat(dx,Nx,1);
h_cell(1) = dx/2;
h_cell(end) = dx/2;
h_osc=h_cell(osc_inds);


theta_dot_0 = (2*pi)./T;

function dydt = inner(t,y)

theta=y(1:N);
I=y((N+1):end);

flux_r=[-D*diff(I)/dx; pr*I(end)];
% dIdx_r(1) = dIdx_r(1)*2;

flux_l=[pl*I(1); flux_r(1:end-1)];

div_flux = (flux_r-flux_l)/dx;

%end-point cells are of length dx/2
div_flux(1) = div_flux(1)*2;
div_flux(end) = div_flux(end)*2;

dydt = [theta_dot_0 - eps*sense(theta).*I(osc_inds); -p*I - div_flux];
dydt(osc_inds_abs) = dydt(osc_inds_abs) + secrete(theta)./h_osc;
% dydt(osc_inds_abs) = secrete(theta);
end


h=@inner;
end