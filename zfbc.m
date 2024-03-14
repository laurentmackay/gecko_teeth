function h=zfbc(ndim)
z =  zeros(ndim,1);
o = ones(ndim,1) ;
function [pl,ql,pr,qr] = inner(xl,ul,xr,ur,t)
pl = z;
pr = z;
ql = o; 
qr = o; 
end
h=@inner;
end