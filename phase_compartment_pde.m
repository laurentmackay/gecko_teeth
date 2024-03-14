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
end



beta=100;

% Nx=N+(N+1)*N_inter;


[xmesh, IJUMP] = oscillator_grid(N, N_inter, L, pad=opts.pad, center='interface');
Nx=length(xmesh);

dx=L/(Nx-1);
T=T0*(1-grad*(0:(N-1))'/(N-1));
% IJUMP_abs = IJUMP+N;


function dydt = inner(t,SOL)

theta = SOL(1:N);
v = SOL(N+1:2*N);
Inh = SOL(2*N+1:end);
Iave = (Inh(IJUMP)+Inh(IJUMP+1))/2;

compflux = beta*(Iave-v);

Icomp =  Iave - dx*compflux/(4*D);    



f_osc = 2*pi./T - eps*sense(theta).*v;
% g_osc =  secrete_prime(theta)*2*pi./T ;
g_osc =  secrete(theta) ;

flux_r=-D*[Inh(2:Nx)-Inh(1:Nx-1); 0]/dx;
flux_l=[0; flux_r(1:end-1)];

flux_r(IJUMP) = -2*D*(Icomp-Inh(IJUMP))/dx;
flux_l(IJUMP+1) = -2*D*(Inh(IJUMP+1)-Icomp)/dx;


 % 0 = -(J(x^+) - J(x^-)) - beta*(v(x)-v)  ::  Flux Jump Condition


div_flux = (flux_r-flux_l)/dx;
% div_flux(1) = div_flux(1)*2;
% div_flux(end) = div_flux(end)*2;



dydt = [f_osc; g_osc + compflux;  -p*Inh - div_flux];



end


h=@inner;
end