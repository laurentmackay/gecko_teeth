function h = phase_compartment_pde( eps, T0,  grad,  D, p, L, N, N_inter, opts)
arguments
    eps  (1,1) {mustBeNumeric} 
    T0  (1,1) {mustBeNumeric} 
    grad  (1,1) {mustBeNumeric} 
    D  (1,1) {mustBeNumeric} 
    p (1,1) {mustBeNumeric} 
    L  (1,1) {mustBeNumeric} 
    N  (1,1) {mustBeInteger} 
    N_inter  (1,1) {mustBeInteger}
    opts.BC (1,:) string = 'Neumann'
    opts.pad (1,1) {mustBeInteger} = N_inter/2
    opts.grad_dist_func (1,1) function_handle = @(x) x
    opts.ablate (:,1) {mustBeNumericOrLogical} = false(N,1)
end

run = ~opts.ablate;


% Nx=N+(N+1)*N_inter;


[xmesh, IJUMP] = oscillator_grid(N, N_inter, L, pad=opts.pad, center='interface');
Nx=length(xmesh);

dx=L/Nx;
x=(0:(N-1))'/(N-1);
T=T0*(1+grad*(opts.grad_dist_func(x-0.5)));

omega_0 =  2*pi./T;

function dydt = inner(t,SOL)

x = SOL(1:N);
y = SOL(N+1:2*N);

r=sqrt(x.^2+y.^2);
theta=atan2(y,x);

% v = SOL(2*N+1:3*N);
Inh = SOL(2*N+1:end);
Iave = (Inh(IJUMP)+Inh(IJUMP+1))/2;

compflux = -run.*secrete(theta);% beta*(Iave-v);

Icomp =  Iave - dx*compflux/(4*D);    



theta_dot = run.*(omega_0 - eps*sense(theta).*Iave);
root = (1-r.^2);


flux_r=-D*[Inh(2:Nx)-Inh(1:Nx-1); 0]/dx;
flux_l=[0; flux_r(1:end-1)];

flux_r(IJUMP) = -2*D*(Icomp-Inh(IJUMP))/dx;
flux_l(IJUMP+1) = -2*D*(Inh(IJUMP+1)-Icomp)/dx;


 % 0 = -(J(x^+) - J(x^-)) - beta*(v(x)-v)  ::  Flux Jump Condition


div_flux = (flux_r-flux_l)/dx;
% div_flux(1) = div_flux(1)*2;
% div_flux(end) = div_flux(end)*2;



dydt = [root.*x - y.*theta_dot; root.*y + x.*theta_dot;  -p*Inh - div_flux];



end


h=@inner;
end