function IC = burnin(rhs,IC, opts)
arguments
    rhs function_handle
    IC  (:,:) {mustBeNumeric}
    opts.tburn (1,1){mustBeNumeric}= 1e4;
    opts.tstep  (1,1) =1e4;
    opts.options struct = odeset()
end
 tburn = opts.tburn;
 tstep = opts.tstep;
%%
tstep=min([tstep, tburn]);
for i=1:tburn/tstep
    [~,sol] = ode15s(rhs , [0,tstep], IC, opts.options );
    IC=sol(end,:);
end
end