
% xosc=xosc(2:end-1);

N0=4;
N=2*(N0)+1;
mid=N0+1;



% T=30;


eps=3.5;
gamma=0.1;
m=1;

x=(0:N-1)/(N-1);

% T=T*(1+0*2*abs(x-0.5) )';

IC=randn(N*2,1)*0.1*1;
odd = [1:2:N, N+1:2:2*N];
even = [2:2:N-1, N+2:2:2*N-1];
IC(even(:))=sqrt(2)/2;
IC(odd(:))=-sqrt(2)/2;

% theta0=randn(N,1);
% 



theta0=zeros(N,1)+0.0*randn(N,1);
% theta0(2:2:end-1)=theta0(2:2:end-1)+pi;

chi=0.95;
adv=0.0;
theta_r = phase_asym_IC3(N0, -0.0);
% theta_r = theta_r - 0.05*randn(1,N0);
theta0(N0+1) = 0;
theta0(N0+2:end) =  theta_r;
theta0(1:N0) =  fliplr(theta_r);

% theta0 = theta0 - 0.05*randn(N,1);

% theta0(1:2:end)=theta0(1:2:end)+pi;

% theta0=theta0+0.01*randn(size(theta0));
% 
% thetal= mod(atan2(ICl(N+1:end),ICl(:,1:N)),2*pi);
% thetar= mod(atan2(ICr(N+1:end),ICr(:,1:N)),2*pi);
% 
% theta0(1:mid-1)= thetal(1:mid-1);
% theta0(mid+1:end)=thetar(mid+1:end);
% theta0(mid)=mod(pi+(thetar(mid+1)+thetal(mid-1))/2,2*pi);

IC(1:N)=cos(theta0);
IC(N+1:end)=sin(theta0);


Jpattern=[speye(N) speye(N); speye(N) speye(N)];

options=odeset('AbsTol',1e-9,'RelTol',1e-9, 'MaxStep', 2*pi/4, 'Jpattern', Jpattern);

rhs = cartesian_phase_inhibitor_ode_nondim(N, gamma, eps);
rhs_back = cartesian_phase_inhibitor_ode_nondim(N, gamma, 0, T0=-1);

%%
tic;
tburn=0;
tburn=1e4;
tstep=2e3;
% tburn=14e3;
% tstep=5e2;
ICP=IC;

% IC=IC_acc';
for i = 1:tburn/tstep
[~,sol] = ode15s(rhs, [0, tstep], ICP , options);
ICP=sol(end,:);
end
if tburn<tstep
    sol=IC';
end
%%
back=false
tmesh=0:0.01:180;
if eps==0 && back
    [t_back,sol_back] = ode15s(@(t,y) rhs_back(t,y), tmesh, sol(end,:) , options);
elseif eps~=0 && back
        [~,sol] = ode15s(@(t,y) rhs_back(t,y), [0,60], sol(end,:) , options);
end
[t_for,sol_for] = ode15s(rhs, tmesh, sol(end,:) , options);

toc 
if eps==0 && back
    tmesh = [t_back(1:end-1); t_for+t_back(end)];
    sol=[flipud(sol_back(2:end,:)); sol_for];
else
    sol=sol_for;
end
% sol=flipud(sol_back);

%%
figure(1);
tskip=0.0*tmesh(end);
iskip = find(tmesh>tskip,1);

sol_theta = mod(atan2(sol(:,N+1:end),sol(:,1:N)),2*pi);
theta = atan2(sol(:,N+1:end),sol(:,1:N));
if back
    tburn =-t_for(end)+2*pi;
end
tburn=0;
inds_full=max(sol_theta)-min(sol_theta)>6
tmin=16;
tburn=-tmin;
T_mid = oscillations_summary(sol_theta, tmesh, tmin=tmin, tmax=1*tmesh(end), dotSize=25, ...
    MinPeakDistance=0, chi_thresh=0.00,  MinPeakProminence=0.95, MinPeakHeight=6, border=true,  cycle_thresh=1.75, ...
    xlines=[0], phi_scatter=[], tburn=tburn, edmund_lines=true, checker=false, period=true, pad=0.02, raster=false, title=false);
% display_oscillations(sol_theta(:,:), tmesh, 0,  MinPeakProminence=0.75)
% yline(0)
% set(gca, 'YTick',[0,1,2,3])
if eps==0
hold on
n=-N0:N0;
T0=1;
dr=mod([sol_theta(1,1:end-1)-sol_theta(1,2:end), nan],2*pi)/(2*pi);
dl=mod([nan,sol_theta(1,2:end)-sol_theta(1,1:end-1)],2*pi)/(2*pi);

% dr0=mod(-[theta0(1:end-1)'-theta0(2:end)', nan],2*pi)/(2*pi);
% dl0=mod(-[nan,theta0(2:end)'-theta0(1:end-1)'],2*pi)/(2*pi);

% if back
%     dr=-[theta0(1:end-1)'-theta0(2:end)', nan]/(2*pi);
%     dl=-[nan,theta0(2:end)'-theta0(1:end-1)']/(2*pi);
% end
yrange=diff(ylim);
pos=get(gca,'Position');
yl=ylim();
for k=0:11

  %   sw =[(1/2).*gamma.^(-1).*k.*((-1)+N).^(-1).*T0.^(-1).*(2.*gamma.*n+(( ...
  % -1)+N).*T0).*(n+(-1).*abs((-1)+n)).^(-1).*(((-1)+N).*T0+2.*gamma.* ...
  % abs((-1)+n))];
  %   plot(n, sw,'-k')

   
     kk=0;
    swc=[(dl-dr+kk*sign(n)).*(T0./(N0*gamma*2*sign(n))).*(N0+gamma*abs(n-1)).*(N0+gamma*abs(1+n))]*0;%+floor(tburn/T_mid)- mod(tburn/T_mid,1);
    h=plot(n, swc,':', 'LineWidth', 2, 'Color',[0,0,0,0.9]);
    angle=180*atan2(diff(pos(4)*swc(end-3:end-2))/yrange,pos(3)/N)/pi;
    text(n(end)-0.5, swc(end-1)+(0.02)*yrange, ['k=' num2str(kk)],'HorizontalAlignment','right','VerticalAlignment','bottom','BackgroundColor',[1,1,1, 0.85],'Clipping','on','Margin',0.5,'Rotation',angle,'FontSize',11)
    uistack(h,"bottom")
    % u

  % for kk=(2*k):(2*k+1)
  % 
  %   swc=[(dl-dr+kk*sign(n)).*(T0./(N0*gamma*2*sign(n))).*(N0+gamma*abs(n-1)).*(N0+gamma*abs(1+n))]*2*pi/T_mid+tburn/T_mid- mod(tburn/T_mid,1);%+floor(tburn/T_mid)- mod(tburn/T_mid,1);
  %   h=plot(n, swc,':', 'LineWidth', 2, 'Color',[0,0,0,0.9]);
  %   angle=180*atan2(diff(pos(4)*swc(end-3:end-2))/yrange,pos(3)/N)/pi;
  %   text(n(end)-0.5, swc(end-1)+(0.02)*yrange, ['k=' num2str(kk)],'HorizontalAlignment','right','VerticalAlignment','bottom','BackgroundColor',[1,1,1, 0.85],'Clipping','on','Margin',0.5,'Rotation',angle,'FontSize',11)
  %   % uistack(h,"bottom")
  % end

% 
% swl=((dl-1/2).*sign(n-1/2)+(k+1/2)).*(T0/(N0*gamma)).*(N0+gamma.*abs(n-1)).*(N0+gamma.*abs(n));
% swl=((dl).*sign(n-1/2)+(k)).*(T0/(N0*gamma)).*(N0+gamma.*abs(n-1)).*(N0+gamma.*abs(n));
% 
% plot(n, swl/T0,'-r')


%   swr=((k+1/2)-(dr-1/2).*sign(n+1/2)).*(T0/(N0*gamma)).*(N0+gamma.*abs(n)).*(N0+gamma.*abs(1+n));
%   swr=((k)-(dr).*sign(n+1/2)).*(T0/(N0*gamma)).*(N0+gamma.*abs(n)).*(N0+gamma.*abs(1+n));
% plot(n, swr,'-b')
% 


% plot(n, abs(((1+k)*((N-1)*T0 + 2*abs(n)*gamma).*((N-1)*T0 +2*gamma+ 2*(abs(n)-1)*gamma))./((4*(N-1)*T0*gamma))),'-r')
end
hold off
end

figure(3); plot(tmesh,abs(cartesian_to_kuramoto(sol)))

return


%%
figure(11);plot(tmesh(tmesh<200), sol_theta(tmesh<200,mid-2),'LineWidth',2)
%%
hr = rectangle('Position',[-N0 190 2*N0 20],'LineStyle','--','EdgeColor',"r");
%%
peaks = get_sol_peaks(sol_theta, tmesh);
[ovec,tvec , phi_l, phi_r, phi_ls, phi_rs] = get_phases(peaks, tmesh);
[chi, chis] = phase_asymmetry(phi_l, phi_r, phi_ls, phi_rs);
[psi, ~] = phase_checker(phi_l, phi_r, phi_ls, phi_rs);
tsample=0:10:3.1e3;
win=15
psi_bar = arrayfun(@(t) mean(psi(tvec>=t-(win/2) & tvec<=t+(win/2)),2,'omitnan'),tsample);
% chi_bar = arrayfun(@(t) mean(chi(tvec>=t-(win/2) & tvec<=t+(win/2)),2,'omitnan'),tsample);
delta_chi_bar = arrayfun(@(t) mean(chi(tvec>=t-(win/2) & tvec<=t+(win/2) & ovec>N0),2,'omitnan')-mean(chi(tvec>=t-(win/2) & tvec<=t+(win/2) & ovec<N0),2,'omitnan'),tsample);

%%
figure(34);
clf();
colororder({'k','r'})

yyaxis left
h1=plot(tsample/(2*pi), psi_bar,'LineWidth',2);
ylabel('$$\bar{\psi}$$','rotation',0,'Interpreter','latex','FontSize',14)
yline(0);
ylim([-0.8, 1])
yyaxis right
h2=plot(tsample/(2*pi), delta_chi_bar,'LineWidth',2);

ylabel('$$\Delta \bar{\chi}$$','rotation',0,'Interpreter','latex','FontSize',14)
hold on;
h3=plot([198,214],[0 0], 'LineWidth',8,'Color','c',LineStyle='-');
hold off;

ylim([-0.8, 1])
hleg = legend([h1,h2,h3],{'anti-phase metric','directionality metric','gecko-like waves'},'FontSize',12, 'location', 'southwest');
% pos = get(hleg,'Position')
% set(hleg,'Position',pos+[0 -1 0 1 ]*0.05)
xlabel('Cycle #');
xlim tight;
%%
figure(33)
plot(t(iskip:end)/(2*pi),sol_theta(iskip:end,1:2), 'LineWidth',2)
legend('\theta_{\rm odd}','\theta_{\rm even}','FontSize',14)
xlabel('t','Interpreter','tex','FontSize',16);
ylabel('\theta','Interpreter','tex','FontSize',16);
set(gca,"YTick",[0,pi, 2*pi], "YTickLabel",["0", "\pi", "2\pi"]);
ylim([0,2*pi])
xlim tight
%exportgraphics(gcf,'anti_phase_ode.pdf')

%%

figure(4)
plot(t(iskip:end)/(2*pi),sol_theta(iskip:end,:), 'LineWidth',2)
legend('\theta_{1}','\theta_{2}','\theta_{3}')
xlabel('t','Interpreter','tex');
ylabel('\theta','Interpreter','tex');
axis tight
%exportgraphics(gcf,'spatially_patterned_ode.pdf')

%%

peaks = get_sol_peaks(sol(:,1:N));


[ovec,tvec,phi_l, phi_r, phi_ls, phi_rs] = get_phases(peaks, tmesh);
[chi, chis] = phase_asymmetry(phi_l, phi_r, phi_ls, phi_rs);
asymmetry_mean = cellfun(@(p) mean(p), chis);
asymmetry_err = cellfun(@(p) std(p), chis);

[psi, psis] = phase_checker(phi_l, phi_r, phi_ls, phi_rs);
checker_mean = cellfun(@(p) mean(p), psis);
checker_err = cellfun(@(p) std(p), psis);
i_left = ovec<N/2;
i_right = ovec>N/2;
figure(8);
scatter(chi(i_left),psi(i_left),5,'filled', 'Color', 'b' );
hold on
scatter(chi(i_right),psi(i_right),5,'filled', 'Color', 'r' )
hold off
ccl=corrcoef(chi(i_left), psi(i_left));
ccl(2)
ccr=corrcoef(chi(i_right), psi(i_right));
ccr(2)  
%%
figure(3)
estimate_periodicity(theta, tmesh, tskip)
%%
disp(u_str(sol(end,:)))
%%
get_per_sol(sol(iskip:end,even(:)), tmesh, sol(end,:), rhs , NTST=600, options=options, periodicity=1, MinPeakProminence=0.75);
%%
fn=['auto/ODE/phase/' num2str(N) '_gamma=' num2str(gamma) '_eps=' num2str(eps) '_antiphase.dat']
write_per_sol(sol(iskip:end,even(:)),tmesh, sol(end,:), rhs , fn, NTST=600, options=options,  periodicity=1, MinPeakProminence=0.75);


%%
 write_ss_sol(sol,'three_phase_anti_locked.dat');

