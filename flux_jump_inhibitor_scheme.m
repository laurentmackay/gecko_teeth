function h = flux_jump_inhibitor_scheme(ode_rhs_and_flux_jumps, kappa, p,  NJUMP, N_inter, opts)
arguments
    ode_rhs_and_flux_jumps  (1,1) function_handle 
    kappa  (1,1) {mustBeNumeric} 
    p (1,1) {mustBeNumeric} 
    NJUMP  (1,1) {mustBeInteger} 
    N_inter  (1,1) {mustBeInteger}
    opts.BC (1,:) string = 'Neumann' %we actually do not handle anything but Neumann right now
    opts.pad (1,1) {mustBeInteger} = N_inter/2
    opts.grad_dist_func (1,1) function_handle = @(x) x
    opts.ablate (:,1) {mustBeNumericOrLogical} = false(NJUMP,1)
    opts.NODE (1,1) {mustBeInteger} = 2
end

NODE = opts.NODE;

alive = ~opts.ablate;
dead = opts.ablate;


% Nx=N+(N+1)*N_inter;

L=NJUMP;

[xmesh, IJUMP] = oscillator_grid(NJUMP, N_inter, L, pad=opts.pad, center='interface');
Nx=length(xmesh);

dx=L/(Nx-1);

kappa2=kappa*kappa;
kappa2dx=kappa2/dx;
D=p*kappa2;
% Ddx = D/dx;
Ddx2 = D/(dx*dx);
rect = [NJUMP,NODE];
i_max_ode = NODE*NJUMP;
vec_sz = [i_max_ode, 1];
GAMMA = zeros(NJUMP,1);
% zero=0;
function dydt = inner(t,SOL)
% zero = zeros(1,size(SOL,2));



ode = reshape(SOL(1:i_max_ode,:), rect);

Inh = SOL(i_max_ode+1:end,:);
Iave = (Inh(IJUMP,:)+Inh(IJUMP+1,:))/2;

[ode, GAMMA] = ode_rhs_and_flux_jumps(ode, Iave, GAMMA);

ode = reshape(alive.*ode, vec_sz);

GAMMA(dead) = 0;

Iint =  Iave + GAMMA/(4*kappa2dx); %inhibitor at the interface, use the ghost point/cell method


flux_r=[-Ddx2*(Inh(2:Nx,:)-Inh(1:Nx-1,:)); 0];
flux_l=[0; flux_r(1:end-1,:)];

flux_r(IJUMP,:) = -2*Ddx2*(Iint-Inh(IJUMP,:));
flux_l(IJUMP+1,:) = -2*Ddx2*(Inh(IJUMP+1,:)-Iint);


dydt = [ode;  -p*Inh - (flux_r-flux_l)];



end


h=@inner;
end