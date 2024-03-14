function v = FHN_u_nullcline(a, opts)
arguments
a {mustBeNumeric}
opts.f0 {mustBeNumeric} = 0
end
v = @(u) -(-1 + u).*u.*(u-a) + opts.f0;
end