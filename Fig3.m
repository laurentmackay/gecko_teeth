%% panel a - run simulation

N0=20;
N=2*(N0)+1;

eps=0.5;
gamma=0;

theta0 = bilateral_phase_waves(N, 0.0);


IC=zeros(N*2,1);
IC(1:N)=cos(theta0);
IC(N+1:end)=sin(theta0);




options=odeset('AbsTol',1e-9,'RelTol',1e-9, 'MaxStep', 2*pi/4, 'Jpattern', [speye(N) speye(N); speye(N) speye(N)]);



rhs = cartesian_phase_inhibitor_ode_nondim(N, gamma, eps);

ICP=burnin(rhs, IC, options=options, tburn=1e4, tstep=2e3);

tmesh=0:0.001:160;

[~,sol] = ode15s(rhs, tmesh, ICP , options);

%% panel a - plotting
figure(1);

sol_theta = mod(atan2(sol(:,N+1:end),sol(:,1:N)),2*pi);
T_mid = oscillations_summary(sol_theta, tmesh, dotSize=25, ...
    MinPeakDistance=0, MinPeakProminence=0.95, MinPeakHeight=6, border=true,  cycle_thresh=1.75, ...
    xlines=0, edmund_lines=false, period=false, pad=0.02, title=false, cycle_max=3.5, border_color= [1,1,1]*0, border_width=0.75);

% exportgraphics(gcf, ['eps=' num2str(eps) '_gamma=' num2str(gamma) '_antiphase_N=' num2str(N)  '.pdf']);


%% panel b - run simulation

N0=20;
N=2*(N0)+1;

eps=8;
gamma=0;

theta0 = bilateral_phase_waves(N, 0.01);
theta0=theta0+randn(size(theta0))

IC=zeros(N*2,1);
IC(1:N)=cos(theta0);
IC(N+1:end)=sin(theta0);



options=odeset('AbsTol',1e-9,'RelTol',1e-9, 'MaxStep', 2*pi/4, 'Jpattern', [speye(N) speye(N); speye(N) speye(N)]);



rhs = cartesian_phase_inhibitor_ode_nondim(N, gamma, eps);

ICP=burnin(rhs, IC, options=options, tburn=1e4, tstep=2e3);

tmesh=0:0.001:160;

[~,sol] = ode15s(rhs, tmesh, ICP , options);

%% panel b - plotting
figure(2);

sol_theta = mod(atan2(sol(:,N+1:end),sol(:,1:N)),2*pi);
T_mid = oscillations_summary(sol_theta, tmesh, dotSize=25, ...
    MinPeakDistance=0, MinPeakProminence=0.95, MinPeakHeight=6, border=true,  cycle_thresh=1.75, ...
    xlines=0, edmund_lines=false, period=false, pad=0.02, title=false, cycle_max=3.5, border_color= [1,1,1]*0, border_width=0.75);

% exportgraphics(gcf, ['eps=' num2str(eps) '_gamma=' num2str(gamma) '_patterned_N=' num2str(N)  '.pdf']);

%% panel c - run simulation

N0=20;
N=2*(N0)+1;

eps=5.75;
gamma=0.07;

theta0 = phase_waves(N, -0.1);
% theta0=theta0+randn(size(theta0))

IC=zeros(N*2,1);
IC(1:N)=cos(theta0);
IC(N+1:end)=sin(theta0);




options=odeset('AbsTol',1e-9,'RelTol',1e-9, 'MaxStep', 2*pi/4, 'Jpattern', [speye(N) speye(N); speye(N) speye(N)]);



rhs = cartesian_phase_inhibitor_ode_nondim(N, gamma, eps);

ICP=burnin(rhs, IC, options=options, tburn=6e4, tstep=2e3);

tmesh=0:0.001:160;

[~,sol] = ode15s(rhs, tmesh, ICP , options);

%% panel c - plotting
figure(3);

sol_theta = mod(atan2(sol(:,N+1:end),sol(:,1:N)),2*pi);
T_mid = oscillations_summary(sol_theta, tmesh, dotSize=35, ...
    MinPeakDistance=0, MinPeakProminence=0.95, MinPeakHeight=6, border=true,  cycle_thresh=1.75, ...
    xlines=0, edmund_lines=true, period=false, pad=0.02, title=false, cycle_max=5, border_color= [1,1,1]*0, border_width=0.75);

% exportgraphics(gcf, ['eps=' num2str(eps) '_gamma=' num2str(gamma) '_rightward_N=' num2str(N)  '.pdf']);


%% panel d - run simulation

N0=45;
N=2*(N0)+1;

eps=6.5;
gamma=0.0;

theta0 = bilateral_phase_waves(N, +0.01);

IC=zeros(N*2,1);
IC(1:N)=cos(theta0);
IC(N+1:end)=sin(theta0);




options=odeset('AbsTol',1e-9,'RelTol',1e-9, 'MaxStep', 2*pi/4, 'Jpattern', [speye(N) speye(N); speye(N) speye(N)]);



rhs = cartesian_phase_inhibitor_ode_nondim(N, gamma, eps);

ICP=burnin(rhs, IC, options=options, tburn=2e4, tstep=2e3);

tmesh=0:0.001:160;

[~,sol] = ode15s(rhs, tmesh, ICP , options);

%% panel d - plotting
figure(4);

sol_theta = mod(atan2(sol(:,N+1:end),sol(:,1:N)),2*pi);

T_mid = oscillations_summary(sol_theta, tmesh, dotSize=25, ...
    MinPeakDistance=0, MinPeakProminence=0.95, MinPeakHeight=6, border=true,  cycle_thresh=1.75, ...
    xlines=0, edmund_lines=true, tburn=-12, period=false, pad=0.02, title=false, border_color= [1,1,1]*0.4, border_width=0.1, cycle_min=0, cycle_max=1.5);

% exportgraphics(gcf, ['eps=' num2str(eps) '_gamma=' num2str(gamma) '_wavelets_N=' num2str(N)  '.pdf']);
