function Jpattern = osc_field_Jpattern(N,N_inter,opt)
arguments
    N (1,1) {mustBeNumeric}
    N_inter (1,1) {mustBeNumeric}
    opt.pad  (1,1) {mustBeInteger} = N_inter/2
    opt.node (1,1) {mustBeInteger} = 1
    opt.npde (1,1) {mustBeInteger} = 1
    opt.center (1,:) string = 'cell'
end

[xmesh, osc_inds] = oscillator_grid(N, N_inter, 1, pad=opt.pad, center=opt.center);

Nx = length(xmesh);
n = opt.node;
m = opt.npde;
J_oo = repmat(speye(N,N),[n,n]);
switch opt.center
    case 'cell'
        J_of = repmat(sparse(1:N,osc_inds', ones(size(osc_inds)), N, Nx), [n,m]);
        J_fo = J_of';
    case 'interface'
        J_of = repmat(sparse([1:N, 1:N],[osc_inds,osc_inds+1 ], ones(size(osc_inds).*[1,2]), N, Nx), [n,m]);
        J_fo = repmat(sparse(1:N,osc_inds', ones(size(osc_inds)), N, Nx), [n,m])';
end
J_ff = repmat(spdiags(ones(Nx,3),-1:1,Nx,Nx),[m,m]);

Jpattern = [J_oo J_of; J_of' J_ff];


end