function h = schic(a,b, eta)
ss=[a ; b/a];
function [u0] = inner(x)
u0=ss.*(1+randn(2,1)*eta);
end
h=@inner;
end