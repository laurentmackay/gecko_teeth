function h = FHN_secrete(a,b, e, k, beta, opts)
arguments
    a  (1,1) {mustBeNumeric} 
    b  (1,1) {mustBeNumeric} 
    e  (:,1) {mustBeNumeric} 
    k  (:,1) {mustBeNumeric} 
    beta  (1,1) {mustBeNumeric}
    opts.f0 (1,1) {mustBeNumeric} = 0;
    opts.timescale (:,1) {mustBeNumeric} = 1;
end



f0=opts.f0;

Tinv = 1./opts.timescale;



function [y, GAMMA] = inner(y, Inh, GAMMA)

u = y(:,1);
v = y(:,2);


GAMMA(:) = beta.*(v-Inh);





f = (u.*(1-u).*(u-a) - v + f0).*Tinv;
g =  e*(k*u-v-b).*Tinv - GAMMA;






%this is a trick that allows for some speed-up apparently
y(:,1) = f;
y(:,2) = g;



end


h=@inner;
end