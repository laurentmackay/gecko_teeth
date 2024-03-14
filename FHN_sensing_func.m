function f = FHN_sensing_func(i, opts)
arguments 
    i (1,1) {mustBeInteger}
    opts.f0 (1,1) {mustBeNumeric} = 0

end
step = smooth_heav();

[x0,n] = FHN_sensing_vectors(abs(i),f0=opts.f0);

if i>0
    f = @(u,v) 1-step(([u,v]-x0)*n,0.0,0.01);
else 
    f= @(u,v) step(([u,v]-x0)*n,0.0,0.01);
end

end