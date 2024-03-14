
load('LG190Raw.mat')
mid=44;
data = cleanupLG190(data, mid);
% peak_inds = get_peak_inds(sol_theta(iskip:end,:));
% locations = -N0:N0;
% EruptionTimes = arrayfun(@(i) [repmat(locations(i), [ length(peak_inds{i}),1]), peak_inds{i}], 1:length(peak_inds), 'UniformOutput', false);
% EruptionTimes =  vertcat(EruptionTimes{:});
% EruptionTimes = sortrows(EruptionTimes,2);


inds = min(data(:,1)):max(data(:,1));

peaks = arrayfun(@(i) data(data(:,1)==i,2)' ,inds ,'UniformOutput' ,false);


%%


% legend({'Simulation','','Data'})

%%
% mid=ceil(length(peaks)/2);
mid=44
good=true(size(inds));
good((mid+8):(mid+22))=false;
% good(mid+11)=true;
% good((mid+11))=true;
good(72-7)=true;
tburn=0;
figure(10);
dim=[700 400];
oscillations_summary([], [], peaks=peaks,edmund_lines=true, period=false, tmin=tburn, tburn=tburn, tmax=250, xlines=[0], border=true, chi_thresh=0.00, checker=false, edmund_flip=false, dotSize=25, cycle_thresh=0.75, xlabel='Tooth #', dim = dim);

[ovec,tvec,phi_l, phi_r, phi_ls, phi_rs] = get_phases(peaks);

bad_rect=[inds(find(~good,1))-mid-2 .5 find(~good,1,'last')-mid-inds(1)/2 3];
rectangle('Position',bad_rect, 'LineStyle','--')

exportgraphics(gcf, 'LG190_asymm.pdf')

figure(11);
oscillations_summary([], [], peaks=peaks,edmund_lines=true, period=false, tmin=tburn, tburn=tburn, tmax=250, xlines=[0], border=true, chi_thresh=0.00, checker=true, edmund_flip=false, dotSize=25, cycle_thresh=0.75, xlabel='Tooth #',  dim = dim);

[ovec,tvec,phi_l, phi_r, phi_ls, phi_rs] = get_phases(peaks);
rectangle('Position',bad_rect, 'LineStyle','--')

exportgraphics(gcf, 'LG190_anti-phase.pdf')

%%


% good(70-8)=true;
% good([67:69]-7)=true;
% good(68-8)=true;

peak_times = cellfun(@(p) p, peaks,'UniformOutput', false);
periods = cellfun(@(t) diff(t), peak_times,'UniformOutput', false);
mean_period = cellfun(@(t) mean(t), periods);
period_err = cellfun(@(t) std(t), periods);

figure(12);
errorbar(inds(good), mean_period(good), period_err(good));
ylabel('Period (days)')
xlabel('Position')

%%
tmesh= 1: max(data(:,2));
periods = cellfun(@(t) diff(tmesh(t)), peaks,'UniformOutput', false);

early = 140;
late = 170;
mid = 44;

mean_period_early = cellfun(@(p,t) mean(p(t(2:end)<=early)),periods, peak_times);
period_err_early = cellfun(@(p,t) std(p(t(2:end)<=early)),periods, peak_times);

mean_period_late = cellfun(@(p,t) mean(p(t(2:end)>=late)),periods, peak_times);
period_err_late= cellfun(@(p,t) std(p(t(2:end)>=late)),periods, peak_times);

figure(13);

ncol=2;
ncol2=1;
subplot(1,ncol,1:(ncol-ncol2))

sym=true;
[ls,gofs]=fit_windowed_periods(peak_times, tmesh, [[0, early];[late, max(data(:,2))]], sym, mask=good);
h1=errorbar(inds(good)-mid, mean_period_early(good), period_err_early(good),'-k','marker','.','LineStyle','none');

hold on

h2=errorbar(inds(good)-mid, mean_period_late(good), period_err_late(good),'-r','marker','.','LineStyle','none');
% h2=errorbar(inds-mid, mean_period_late, period_err_late,'-r','marker','.','LineStyle','none');
l0=ls{1};
l1=ls{2};
plot(inds-mid, l0(inds-inds(1)),'-k','LineWidth',2)



plot(inds-mid, l1(inds-inds(1)),'-r','LineWidth',2)

hold off
ylabel('Period (days)')
xlabel('Tooth #')
set(gca,'FontSize',12)

[~,~,DeltaT,~] = get_pars_and_CIs(ls);
if sym


    gamma_early = DeltaT(1);
    gamma_late = DeltaT(2);


else
    gamma_left_early = l0.gamma_l;
    gamma_right_early = l0.gamma_r;
    gamma_early = mean([gamma_left_early, gamma_right_early]);

    gamma_left_late = l1.gamma_l;
    gamma_right_late = l1.gamma_r;
    gamma_late  = mean([gamma_left_late, gamma_right_late])
end


cycle_map=peaktimes_to_cycles(peak_times(mid),tmesh);

label_early = [ 'Cycle # < ' num2str(round(cycle_map(early),1)) ' (\DeltaT\approx' num2str(gamma_early,3) ')'];
label_late = [ 'Cycle # > ' num2str(round(cycle_map(late),1)) ' (\DeltaT\approx' num2str(gamma_late,3) ')'];

legend([h1, h2],label_early,label_late,'FontSize',12)


% figure(14);
win_size=150;
wins =cell2mat(arrayfun(@(t) [t-win_size/2, t+win_size/2], linspace(1, max(data(:,2)),60) ,'UniformOutput',false)');



% [ls,gofs]=fit_windowed_periods(peak_times, tmesh, wins, sym, mask=good, plot=2);
plotPeriodGradientDynamics(peak_times, tmesh, wins, ax1=subplot(2,ncol,(ncol-ncol2+1):ncol), ax2=subplot(2,ncol,ncol+((ncol-ncol2+1):ncol)), mask=good, plotfits=false)

set(13,'Position',[321   325   800   500])
% tightFig(13)

return

%%
% figure(2);surf(times, inds, mean_period_win');
% 
% hold on
% surf(times, inds, cell2mat(cellfun(@(cfit) cfit(inds),ls,'UniformOutput',false)),'EdgeColor','none','FaceColor',[0,0,0])
% title(['mean r^2=' num2str(mean(cellfun(@(gof) gof.rsquare, gofs))) ])
% hold off