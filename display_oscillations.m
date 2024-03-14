function  display_oscillations(sol, tmesh, tskip, extra)
arguments
    sol (:,:) {mustBeNumeric}
    tmesh (1,:) {mustBeNumeric}
    tskip (1,1) {mustBeNumeric}
    extra.fields (1,:) cell = []
    extra.nrows  (1,1) {mustBeInteger} = 3
    extra.checker (1,1) {mustBeNumericOrLogical} = false;
    extra.MinPeakDistance (1,1) {mustBeNumeric} = 0;
    extra.MinPeakProminence (1,1) {mustBeNumeric} = 0.1;
end

upper_style = {'XTick',[], 'FontSize',14};
lower_style = { 'FontSize',18};

iskip=find(tmesh>=tskip,1);
sol=sol(iskip:end,:);
tmesh=tmesh(iskip:end);

N=size(sol,2);
if N>1
    mid=ceil(N/2)+(1-mod(N,2))/2;
    oscs = (1:N)-mid;
else
    mid=1;
    oscs=[0];
end

base_panels = 2+extra.checker;
N_panels=base_panels+length(extra.fields);
nrows = extra.nrows;
colspan= [nrows-1,1];
tiledlayout( nrows, N_panels);



peaks = get_sol_peaks(sol, tmesh, MinPeakDistance=extra.MinPeakDistance, MinPeakProminence = extra.MinPeakProminence);

peak_times = cellfun(@(p) p, peaks,'UniformOutput', false);
mean_period = cellfun(@(t) mean(diff(t)), peak_times);
period_err = cellfun(@(t) std(diff(t)), peak_times);


nexttile(1);
errorbar(oscs, mean_period, period_err);
ylabel('Period (days)', 'FontSize',14)

set(gca, upper_style{:})


nexttile(N_panels+1,colspan);
colormap(gca, 'parula')
s=imagesc(oscs, flipud(tmesh/365),flipud(sol));
set(gca,'YTickLabel',flipud(get(gca, 'YTickLabel')));
yt0=get(gca,'YTick');
set(gca,'YTick',yt0-(yt0(1)-tmesh(1)/365));
% s.EdgeAlpha=0;




colorbar
% colormap(gca, 'parula');

xlabel('Oscillator #');
ylabel('Time (years)')
set(gca,  lower_style{:})


[ovec,tvec,phi_l, phi_r, phi_ls, phi_rs] = get_phases(peaks);
% tvec=tmesh(tvec);
[chi, chis] = phase_asymmetry(phi_l, phi_r, phi_ls, phi_rs);
asymmetry_mean = cellfun(@(p) mean(p,'all','omitnan'), chis);
asymmetry_err = cellfun(@(p) std(p,[],'omitnan'), chis);

nexttile(2);
errorbar(oscs, asymmetry_mean, asymmetry_err)

yline(0);

ylabel({'Phase',  'Asymmetry'});
set(gca, upper_style{:})

nexttile(N_panels+2, colspan);


scatter(ovec-mid, tvec/365, 55, chi, "filled");
xlabel('Oscillator #');
ylabel('Time (years)')

set(gca, lower_style{:})

xlim([1,N]-mid)
pad=0.01;
range=(tmesh(end)-tmesh(1))/365;
ylim(tmesh(1)/365+[-pad*range (1+pad)*range])

% map=[127,191,123
%     247,247,247;
%     175,141,195]/255;
%
% colormap(gca,padmap(map,25))
map = customcolormap_preset('purple-white-green');
map=sqrt(map);
colormap(gca, map);
colorbar

% clim([-1, 1])
clim([-0.2, 0.2])

if extra.checker
[psi, psis] = phase_checker(phi_l, phi_r, phi_ls, phi_rs);
checker_mean = cellfun(@(p) mean(p), psis);
checker_err = cellfun(@(p) std(p), psis);

nexttile(base_panels);
errorbar(oscs, checker_mean, checker_err)

yline(0);

ylabel('Checkered');
set(gca, upper_style{:})

nexttile(N_panels+base_panels, colspan);


scatter(ovec-mid, tvec/365, 55, psi, "filled");
xlabel('Oscillator #');
ylabel('Time (years)')

set(gca, lower_style{:})

xlim([1,N]-mid)
pad=0.01;
range=(tmesh(end)-tmesh(1))/365;
ylim(tmesh(1)/365+[-pad*range (1+pad)*range])

colormap(gca, hot);
clim([0,1])
colorbar
end


for i=1:length(extra.fields)
    field = extra.fields{i};
    field = field(iskip:end,:);
    x=1:size(field,2);
    nexttile(base_panels+i);

    h=plot(x,[mean(field,1); max(field,[],1); min(field,[],1)], 'Color','k');
    set(h(2:3),'LineStyle',':')
    % h=plot(x,mean(field,1), 'Color','k');
    %

    xlim([1, size(field,2)])


    nexttile(N_panels+base_panels+i, colspan);

    s=imagesc(x, flipud(tmesh), flipud(field));
    set(gca,'YTickLabel',flipud(get(gca, 'YTickLabel')));
    yt0=get(gca,'YTick');
    set(gca,'YTick',yt0-(yt0(1)-tmesh(1)));
    colorbar
end

end