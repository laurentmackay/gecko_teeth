function h = FHN_secrete(a,b, e, k, eps, opts)
arguments
    a  (1,1) {mustBeNumeric} 
    b  (1,1) {mustBeNumeric} 
    e  (:,1) {mustBeNumeric} 
    k  (:,1) {mustBeNumeric} 
    eps (1,1) {mustBeNumeric}
    opts.f0 (1,1) {mustBeNumeric} = 0;
    opts.timescale (:,1) {mustBeNumeric} = 1;
    opts.sense_phase  (1,1) {mustBeInteger}=1;
    opts.secrete_phase  (1,1) {mustBeInteger}=4;
end



f0=opts.f0;

Tinv = 1./opts.timescale;

% m=12;
% step = smooth_heav();

% [x0,n] = FHN_sensing_vectors(opts.sense_phase,f0=opts.f0);
% [x02,n2] = FHN_sensing_vectors(opts.secrete_phase,f0=opts.f0);
% 
% sense = @(u,v) 1-step(([u,v]-x0)*n,0.0,0.01);
% secrete = @(u,v) 1-step(([u,v]-x02)*n2,0.0,0.01);
sense = FHN_sensing_func(opts.sense_phase,f0=opts.f0);
function [y, GAMMA] = inner(y, Inh, GAMMA)

u = y(:,1);
v = y(:,2);


GAMMA(:) = v;





f = (u.*(1-u).*(u-a) - v + f0).*Tinv-eps*sense(u,v).*Inh;
% f = (u.*(1-u).*(u-a) - v + f0).*Tinv-alpha.*Inh;
g =  e*(k*u-v-b).*Tinv;
% f=f.*(1-alpha*sense(u,v).*Inh);
% g=g.*(1-alpha*sense(u,v).*Inh);





%this is a trick that allows for some speed-up apparently
y(:,1) = f;
y(:,2) = g;



end


h=@inner;
end