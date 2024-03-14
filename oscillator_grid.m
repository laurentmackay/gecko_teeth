function [xmesh,osc_inds] = oscillator_grid(N, N_inter, L, opts)
arguments
    N  (1,1) {mustBeInteger} 
    N_inter  (1,1) {mustBeInteger}
    L  (1,1) {mustBeNumeric} 
    opts.pad (1,1) {mustBeInteger} = N_inter
    opts.center (1,:) string = 'cell'
end
pad = opts.pad;




Nx=(N-1)*(N_inter)+2*pad;
switch opts.center

    case 'cell'
        Nx=Nx+N;
        d=(N_inter+1);
        osc_inds = (pad+1:d:(Nx-pad))';
        xmesh=linspace(0,L,Nx)';
        
    case 'interface'
        osc_inds = pad:N_inter:Nx-pad;
        dx=L/Nx;
        xmesh=linspace(dx/2,L-dx/2,Nx)';
    otherwise
        warning('unexpeted center value')
        osc_inds=[];

end


end