function h = schpde(a,b,Du,Dv)

c0=[1;1];
D=[Du; Dv];
function [c,f,s] = inner(x,t,u,dudx)
c=c0;
f=D.*dudx;
rx=u(1)*u(1)*u(2);
s=[(a-u(1)+rx); (b-rx)];
end


h=@inner;
end