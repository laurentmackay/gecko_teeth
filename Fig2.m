
%% panel a - run simulations
N0=21;
N=2*(N0)+1;

eps=0;
gamma=1.0;

theta0 = bilateral_phase_waves(N, 0.0);


IC=zeros(N*2,1);
IC(1:N)=cos(theta0);
IC(N+1:end)=sin(theta0);



options=odeset('AbsTol',1e-9,'RelTol',1e-9, 'MaxStep', 2*pi/4, 'Jpattern', [speye(N) speye(N); speye(N) speye(N)]);

rhs = cartesian_phase_inhibitor_ode_nondim(N, gamma, 0);
rhs_back = cartesian_phase_inhibitor_ode_nondim(N, gamma, 0, T0=-1);

tmesh=0:0.001:60;

[t_back,sol_back] = ode15s(@(t,y) rhs_back(t,y), tmesh,IC , options);
[t_for,sol_for] = ode15s(rhs, tmesh, IC , options);


tmesh = [t_back(1:end-1); t_for+t_back(end)];
sol=[flipud(sol_back(2:end,:)); sol_for];


%% panel a - plotting

figure(1);
clf();
sol_theta = mod(atan2(sol(:,N+1:end),sol(:,1:N)),2*pi);
T_mid = oscillations_summary(sol_theta, tmesh, tmax=1*tmesh(end), dotSize=25, ...
    MinPeakDistance=0, chi_thresh=1e-3,  MinPeakProminence=0.95, MinPeakHeight=6, border=true,  cycle_thresh=1.75, ...
    xlines=0, tburn=-t_back(end), edmund_lines=true, checker=false, period=true, pad=0.02, title=false, ...
    cycle_min=-5, cycle_max=5);


xx=0.03;
y1=0.02;
y2=0.5;
annotation('textbox',[xx y1 .0 .0],'String','back-to-front','LineStyle','none','Rotation',90,'FontSize',15);
annotation('textbox',[xx y2 .0 .0],'String','front-to-back','LineStyle','none','Rotation',90,'FontSize',15);
apos=get(gca, 'Position');
d=0.046;
c=[1,1,1]*0.7;
annotation('line',[pos(1) xx+d],[apos(2) y1+0.015],'LineWidth',0.5, 'Color',c);
annotation('line',[pos(1) xx+d],[apos(2)+apos(4)/2-0.001 y1+0.24],'LineWidth',0.5, 'Color',c);

annotation('line',[pos(1) xx+d],[apos(2)+apos(4) y2+0.24],'LineWidth',0.5, 'Color',c);
annotation('line',[pos(1) xx+d],[apos(2)+apos(4)/2+0.001 y2+0.015],'LineWidth',0.5, 'Color',c);


plot_switching_times(N0, 2*pi, sol_theta(tmesh==0,:), gamma, bottom=true, k=0, tburn=-t_back(end))

% exportgraphics(gcf, ['uncoupled_antiphase_switch_gamma=' num2str(gamma) '.pdf'],'ContentType','vector')



%% panel b - run simulation
% most of this is the same as for panel a

N0=21;
N=2*(N0)+1;

eps=0;
gamma=1.0;

theta0 = bilateral_phase_waves(N, 0.0);


IC=zeros(N*2,1);
IC(1:N)=cos(theta0);
IC(N+1:end)=sin(theta0);



options=odeset('AbsTol',1e-9,'RelTol',1e-9, 'MaxStep', 2*pi/4, 'Jpattern', [speye(N) speye(N); speye(N) speye(N)]);

rhs = cartesian_phase_inhibitor_ode_nondim(N, gamma, 0);

tmesh=0:0.001:600;

[~,sol] = ode15s(rhs, tmesh, IC , options);

%% panel b - plotting

figure(2);
sol_theta = mod(atan2(sol(:,N+1:end),sol(:,1:N)),2*pi);
T_mid = oscillations_summary(sol_theta, tmesh, tmax=1*tmesh(end), dotSize=25, ...
    MinPeakDistance=0, MinPeakProminence=0.95, MinPeakHeight=6, border=false,  cycle_thresh=1.75, ...
    xlines=0, edmund_lines=false, period=true, pad=0.02, title=false);


plot_switching_times(N0, 1, sol_theta(1,:), gamma)


% exportgraphics(gcf, ['uncoupled_longterm_gamma=' num2str(gamma) '.pdf'],'ContentType','vector')


%% panel c - run simulation
% most of this is the same as for panel a

N0=45;
N=2*(N0)+1;

eps=0;
gamma=0.1;

theta0 = bilateral_phase_waves(N, 0.0, delta0=0);


IC=zeros(N*2,1);
IC(1:N)=cos(theta0);
IC(N+1:end)=sin(theta0);



options=odeset('AbsTol',1e-9,'RelTol',1e-9, 'MaxStep', 2*pi/4, 'Jpattern', [speye(N) speye(N); speye(N) speye(N)]);

rhs = cartesian_phase_inhibitor_ode_nondim(N, gamma, 0);

tmesh=0:0.01:6400;

[~,sol] = ode15s(rhs, tmesh, IC , options);

%% panel c - plotting

figure(3);
sol_theta = mod(atan2(sol(:,N+1:end),sol(:,1:N)),2*pi);
T_mid = oscillations_summary(sol_theta, tmesh, dotSize=25, ...
    MinPeakDistance=0, MinPeakProminence=0.95, MinPeakHeight=6, border=false,  cycle_thresh=1.75, ...
    xlines=0, edmund_lines=false, raster=true, period=true, pad=0.02, title=false);


plot_switching_times(N0, 1, sol_theta(1,:), gamma)

% exportgraphics(gcf, ['uncoupled_longterm_gamma=' num2str(gamma) '.pdf'],'ContentType','vector')