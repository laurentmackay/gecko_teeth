function  write_per_sol(samp, t, IC, rhsfun, path, opts)
arguments
    samp  (:,:) {mustBeNumeric} 
    t  (:,1) {mustBeNumeric} 
    IC  (:,1) {mustBeNumeric} 
    rhsfun (1,1) function_handle
    path (1,:) char
    opts.options (1,1) struct = odeset()
    opts.NTST (1,1) {mustBeNumeric}  = 150
    opts.periodicity  (1,1) {mustBeInteger} = 1
    opts.MinPeakProminence (1,1) {mustBeNumeric} = 0.1
end
[sol, T]=get_per_sol(samp, t, IC, rhsfun, NTST=opts.NTST, options=opts.options, periodicity=opts.periodicity, MinPeakProminence=opts.MinPeakProminence);

fid = fopen( path ,'w','n');

fprintf(fid, [ repmat('%2.9E ', 1, length(IC)+1) '\n'], [T(1:end) sol(1:end,:)]');


% fprintf( [ repmat('%2.9E ', 1, length(IC)+1) '\n'], [T(1:end) sol(1:end,:)]')

end