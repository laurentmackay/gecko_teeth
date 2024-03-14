xmesh=linspace(0,10,500);
tmesh=linspace(0,1000,500);

a=8e-4;
b=.505;
Du=0.001; Dv=1;
r=sqrt(Du/Dv)
acr=b^2*r*(b*r+2)
amax = min(b,b^3)
a>acr && a<amax

sol = pdepe(0,schpde(a,b,Du,Dv), schic(a,b,0.1),zfbc(2),xmesh,tmesh);

figure(1);
imagesc(sol(:,:,1));
colorbar