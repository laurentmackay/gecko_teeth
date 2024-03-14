xmesh=linspace(0,160,800);
tmesh=linspace(0,140,350);

a=3;
b=10.445;
sigmac = ((sqrt(1+a^2)-1)/a)^2;
Dv=10;

Du=1.1*sigmac*Dv
r=sqrt(Du/Dv);
bhopf=1+a^2
bturing = (1+a*r)^2
a + (2*r)/(-1 + r^2)


sol = pdepe(0,brpde_grad(a,b,Du,Dv, xmesh(end)), bric(a,b,0.05), zfbc(2), xmesh, tmesh);
% ic = cell2mat(arrayfun(bric(a,b,0.01), xmesh, 'UniformOutput', false));
% sol = pbcpdeSolver(brpde_pbc(a,b,Du,Dv), ic, xmesh, tmesh);
%%
figure(1);
iskip=30;
s=pcolor(xmesh, tmesh(iskip:end),sol(iskip:end,:,2));
s.EdgeAlpha=0;
xlabel('space','FontSize',16);
ylabel('time','FontSize',16);
colorbar