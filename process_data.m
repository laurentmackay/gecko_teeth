
load('LG190Raw.mat')
data = cleanupLG190(data, 44);
% peak_inds = get_peak_inds(sol_theta(iskip:end,:));
% locations = -N0:N0;
% EruptionTimes = arrayfun(@(i) [repmat(locations(i), [ length(peak_inds{i}),1]), peak_inds{i}], 1:length(peak_inds), 'UniformOutput', false);
% EruptionTimes =  vertcat(EruptionTimes{:});
% EruptionTimes = sortrows(EruptionTimes,2);


inds = min(data(:,1)):max(data(:,1));

peaks = arrayfun(@(i) data(data(:,1)==i,2) ,inds ,'UniformOutput' ,false);
figure(1);
tburn=0;
[~,tmap] = oscillations_summary([], [], peaks=peaks,edmund_lines=true, period=false, tmin=tburn, tburn=tburn, tmax=250, xlines=[0], border=true, chi_thresh=0.00, checker=false, edmund_flip=false, dotSize=35, cycle_thresh=0.75);

[ovec,tvec,phi_l, phi_r, phi_ls, phi_rs] = get_phases(peaks);
rectangle('Position',[8 .75 23 3], 'LineStyle','--')

%%



[x,y] = meshgrid(0:0.01:1,0:0.01:1);
chi = (sin(pi*(x-y))/pi);
chi = (x-y);
psi = (sin(pi*x).*sin(pi*y)).^(3);


phases = [cell2mat([phi_rs{:}]); cell2mat([phi_ls{:}])];
r=sqrt(sum((phases-[0.5;0.5]).^2));
rsamp = 0:0.001:sqrt(2)/2;
[pdf, ~] = ksdensity(r, rsamp,'function','cdf','BoundaryCorrection','reflection');
figure(52);
hold on
plot(rsamp,1-(pdf)/max(pdf),rsamp, sin(pi*(rsamp/sqrt(2)+0.5)).^(2))
plot(rsamp,1-(pdf)/max(pdf),rsamp, sin(pi*(rsamp/sqrt(2)+0.5)).^(2*2))
plot(rsamp,1-(pdf)/max(pdf),rsamp, sin(pi*(rsamp/sqrt(2)+0.5)).^(2*3))
hold off
%%
figure(11);
subplot(1,3,1);
s=pcolor(x,y,chi);
set(s,'LineStyle','none')
colorbar
hold on
h=phi_scatter(cell2mat([phi_rs{:}]), cell2mat([phi_ls{:}]),'filled');
hold off







subplot(1,3,2);
s=pcolor(x,y,psi);
set(s,'LineStyle','none')
colorbar
hold on
h=phi_scatter(cell2mat([phi_rs{:}]), cell2mat([phi_ls{:}]),'filled');
hold off

subplot(1,3,3);
s=pcolor(x,y,chi.*psi);
set(s,'LineStyle','none')
colorbar
hold on
h=phi_scatter(cell2mat([phi_rs{:}]), cell2mat([phi_ls{:}]),'filled');
hold off

set(gcf, 'Position', [153         524        1618         420])

% legend({'Simulation','','Data'})
