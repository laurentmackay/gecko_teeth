function h = FHN_oscillator_pde(b, e,  D, p, L, N, N_inter, opts)
arguments
    b  (1,1) {mustBeNumeric} 
    e  (:,1) {mustBeNumeric} 
    D  (1,1) {mustBeNumeric} 
    p (1,1) {mustBeNumeric} 
    L  (1,1) {mustBeNumeric} 
    N  (1,1) {mustBeInteger} 
    N_inter  (1,1) {mustBeInteger}
    opts.BC (1,:) string = 'Neumann'
    opts.pad (1,1) {mustBeInteger} = N_inter/2
end
a = 0.169;
k = 0.6;
% e = 1e-2;

beta=100;

% Nx=N+(N+1)*N_inter;


[xmesh, IJUMP] = oscillator_grid(N, N_inter, L, pad=opts.pad, center='interface');
Nx=length(xmesh);

dx=L/(Nx-1);

% IJUMP_abs = IJUMP+N;


function dydt = inner(t,SOL)

u = SOL(1:N);
v = SOL(N+1:2*N);
Inh = SOL(2*N+1:end);
Iave = (Inh(IJUMP)+Inh(IJUMP+1))/2;
% compflux = beta*(Icomp-v);
compflux = beta*(Iave-v);

Icomp =  Iave - dx*compflux/(4*D);    
f_osc= u.*(1-u).*(u-a)-v;
g_osc =  e.*(k*u-v-b);

flux_r=[-D*(Inh(2:Nx)-Inh(1:Nx-1))/dx; 0];
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