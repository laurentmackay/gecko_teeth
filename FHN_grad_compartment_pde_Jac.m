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
    opts.BC (1,:) string = 'Neumann'
    opts.pad (1,1) {mustBeInteger} = N_inter/2
    opts.grad_dist_func (1,1) function_handle = @(x) x
end
a = 0.169;
k = 0.6;
% e = 1e-2;

beta=100;

% Nx=N+(N+1)*N_inter;


[xmesh, IJUMP] = oscillator_grid(N, N_inter, L, pad=opts.pad, center='interface');
Nx=length(xmesh);

dx=L/(Nx-1);
x=(0:(N-1))'/(N-1);
T=1+grad*(opts.grad_dist_func(x-0.5));
Tinv = 1./T;
% T=1+grad*x;
% IJUMP_abs = IJUMP+N;

Ddx = D/dx;
Ntot=2*N+Nx;

DFUDV = -Tinv';
DFVDU = e*k*Tinv';
DFVDV = -e*Tinv'-beta;
DFVDI = repmat(beta/2, 1, N);
DFIDV = repmat(beta/2*dx, 1, N);

e = ones(Nx,1);
A = spdiags([e -2*e e],-1:1,Nx,Nx);
A(1,1)=-1;
A(Nx,Nx)=1;

RDmat = D*A/(dx^2)-p*speye(Nx,Nx);

for k=1:N
    RDmat(IJUMP(k),IJUMP(k))=RDmat(IJUMP(k),IJUMP(k))-beta/(4*dx);
    RDmat(IJUMP(k),IJUMP(k)+1)=RDmat(IJUMP(k),IJUMP(k)+1)-beta/(4*dx);
    RDmat(IJUMP(k)+1,IJUMP(k))=RDmat(IJUMP(k)+1,IJUMP(k))-beta/(4*dx);
    RDmat(IJUMP(k)+1,IJUMP(k)+1)=RDmat(IJUMP(k)+1,IJUMP(k)+1)-beta/(4*dx);
end

[iRD, jRD, vRD] = find(RDmat);

uind=1:N;
vind = N+uind;
IJUMP_abs = 2*N + IJUMP;
i=    [uind  uind  vind  vind     vind     vind    IJUMP_abs+1 IJUMP_abs+1  2*N+iRD']';
j=    [uind  vind  uind  vind  IJUMP_abs IJUMP_abs+1     vind      vind     2*N+jRD']';
v0 =  [      DFUDV DFVDU DFVDV    DFVDI    DFVDI         DFIDV     DFIDV        vRD']';
DFDY0 = sparse(i,j,[zeros(N,1); v0 ], Ntot, Ntot, length(v0)+N);
% DFDY0 = full(DFDY0);

iuu=sub2ind([Ntot, Ntot], 1:N, 1:N);

function DFDY = inner(t,Y)

u = Y(1:N);



DFUDU = ((1-u).*(u-a) - u.*(u-a) + u.*(1-u)).*Tinv;



% tic

DFDY = sparse(i,j, [DFUDU; v0 ], Ntot, Ntot);
% DFDY=DFDY0;
% DFDY(iuu)=DFUDU;
% disp(t)
% toc




end


h=@inner;
end