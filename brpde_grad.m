function h = brpde_grad(a,b,Du,Dv, L)

c0=[1;1];
D=[Du; Dv];
function [c,f,s] = inner(x,t,u,dudx)
c=c0;
f=D.*dudx;
rx=u(1)*u(1)*u(2);
s=[a-(b+1)*u(1)+rx; (b*u(1)-rx)];
end


h=@inner;
end