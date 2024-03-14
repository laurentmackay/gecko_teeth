function h = schpde(a,b,Du,Dv)

D0=[Du; Dv];
function [D,s] = inner(x,t,u,dudx)

D=D0;
rx=u(1,:).*u(1,:).*u(2,:);
s=[a-(b+1)*u(1,:)+rx; (b*u(1,:)-rx)];
end


h=@inner;
end