function  write_per_sol(sol, name, opts)
arguments
    sol (:,:) {mustBeNumeric}
    name (1,:) char
    opts.options (1,1) struct = odeset()
    opts.NTST (1,1) {mustBeNumeric}  = 150
end


fid = fopen([ '/home/lmackay/atest/' name ],'w','n');

fprintf(fid, [ repmat('%2.9E ', 1, size(sol,2)+1) '\n'], [0.0 sol(end,:)]');


% fprintf( [ repmat('%2.9E ', 1, length(IC)+1) '\n'], [T(1:end) sol(1:end,:)]')

end