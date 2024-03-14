function  plotPeriodGradientDynamics(peaks, tmesh, wins, opts)
arguments

    peaks (1,:) {mustBeA(peaks, 'cell')}
    tmesh (1,:) {mustBeNumeric}
    wins (:,2) {mustBeNumeric}
    opts.sym (1,1) {mustBeNumericOrLogical} = true;
    opts.ax1 (1,1) {mustBeA(opts.ax1, 'matlab.graphics.axis.Axes')} = subplot(1,2,1);
    opts.ax2 (1,1) {mustBeA(opts.ax2, 'matlab.graphics.axis.Axes')} = subplot(1,2,2);
    opts.units (1,:) {mustBeText} = 'days';
    opts.plotfits (1,1) {mustBeNumericOrLogical} = false
    opts.mask (1,:) {mustBeNumericOrLogical} = true(1,length(peaks));
end

[fits,~]=fit_windowed_periods(peaks, tmesh, wins, opts.sym, mask=opts.mask, plot=opts.plotfits);

if opts.sym
    [T0, T0_CI, DeltaT, DeltaT_CI] = get_pars_and_CIs(fits);
else
    DeltaT = cellfun(@(cfit) (cfit.gamma_l+cfit.gamma_r)/2,fits );
    DeltaT_CI = cell2mat(cellfun(@(cfit) (get_col(confint(cfit),2)+get_col(confint(cfit),3))/2,fits ,'UniformOutput',false));
end

ax1=opts.ax1;
ax2=opts.ax2;

mid=ceil(length(peaks)/2);
win_centers = mean(wins,2);
[cycle_map,~]=peaktimes_to_cycles(peaks, tmesh, ind=mid);
cycles=cycle_map(win_centers);



plot(ax1, cycles, DeltaT,'-k', 'LineWidth',2);
hold(ax1,"on");
plot(ax1, cycles, DeltaT_CI,'--k', 'LineWidth',1 )
hold(ax1,"off")
axis(ax1,"tight");
ylim(ax1, [0, max(DeltaT_CI(:))])
set(ax1,'XTick',[])
% xlabel('day','FontSize',14);
ylabel(ax1,['\DeltaT (' opts.units ')'],'FontSize',12)
set(ax1,'FontSize',12)


plot(ax2, cycles,T0,'-k', 'LineWidth',2)
hold(ax2,"on");
plot(ax2, cycles,T0_CI,'--k', 'LineWidth',1 )
hold(ax2,"off");
axis(ax2,"tight");
% ylim([0, max(Tmid_CI(:))])
xlabel(ax2,'Cycle #','FontSize',14);
ylabel(ax2,['T_{0} (' opts.units ')'],'FontSize',12)
set(ax2,'FontSize',12)

end