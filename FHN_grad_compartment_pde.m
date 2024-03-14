function h = FHN_oscillator_pde(b, e, grad,  D, p, L, N, N_inter, opts)
arguments
    b  (1,1) {mustBeNumeric} 
    e  (:,1) {mustBeNumeric} 
    grad  (1,1) {mustBeNumeric} 
    D  (1,1) {mustBeNumeric} 
    p (1,1) {mustBeNumeric} 
    L  (1,1) {mustBeNumeric} 
    N  (1,1) {mustBeInteger} 
    N_inter  (1,1) {mustBeInteger}
    opts.f0 (1,1) {mustBeNumeric} = 0
    opts.BC (1,:) string = 'Neumann'
    opts.pad (1,1) {mustBeInteger} = N_inter/2
    opts.grad_dist_func (1,1) function_handle = @(x) x
    opts.ablate (:,1) {mustBeNumericOrLogical} = false(N,1)
end

f0=opts.f0;

run = ~opts.ablate;

a = 0.169;
k = 0.6;
% e = 1e-2;

beta=0.01;

% Nx=N+(N+1)*N_inter;


[xmesh, IJUMP] = oscillator_grid(N, N_inter, L, pad=opts.pad, center='interface');
Nx=length(xmesh);

dx=L/(Nx-1);
x=(0:(N-1))'/(N-1);
T=1+grad*(opts.grad_dist_func(x-0.5));
Tinv = run./T;
% T=1+grad*x;
% IJUMP_abs = IJUMP+N;

Ddx = D/dx;
function dydt = inner(t,SOL)
zero = zeros(1,size(SOL,2));
% zero=0;
u = SOL(1:N,:);
v = SOL(N+1:2*N,:);
Inh = SOL(2*N+1:end,:);
Iave = (Inh(IJUMP,:)+Inh(IJUMP+1,:))/2;

compflux = beta*run.*(Iave-v);

Icomp =  Iave - compflux/(4*Ddx);    
f_osc= (u.*(1-u).*(u-a)-v +f0 ).*Tinv;
g_osc =  e*(k*u-v-b).*Tinv;

flux_r=[-Ddx*(Inh(2:Nx,:)-Inh(1:Nx-1,:)); zero];
flux_l=[zero; flux_r(1:end-1,:)];

flux_r(IJUMP,:) = -2*Ddx*(Icomp-Inh(IJUMP,:));
flux_l(IJUMP+1,:) = -2*Ddx*(Inh(IJUMP+1,:)-Icomp);






dydt = [f_osc; g_osc + compflux;  -p*Inh - (flux_r-flux_l)/dx];



end


h=@inner;
end