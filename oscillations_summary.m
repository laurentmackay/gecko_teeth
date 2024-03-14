function  T_mid=oscillations_summary(sol, tmesh, opts)
arguments
    sol (:,:) {mustBeNumeric}
    tmesh (1,:) {mustBeNumeric}
    opts.MinPeakDistance (1,1) {mustBeNumeric} = 0;
    opts.MinPeakProminence (1,1) {mustBeNumeric} = 0.1;
    opts.MinPeakHeight (1,1) {mustBeNumeric} = 1e-5;
    opts.border (1,1) {mustBeNumericOrLogical} = false;
    opts.xlines (:,1) {mustBeNumeric} = []
    opts.tmin (1,1) {mustBeNumeric} = -Inf;
    opts.tmax (1,1) {mustBeNumeric} = Inf;
    opts.cycle_min (:,:) {mustBeNumeric} = [];
    opts.cycle_max (:,:) {mustBeNumeric} = [];
    opts.clim (:,1) {mustBeNumeric} = [-0.2,0.2]
    opts.phi_scatter (:,:) {mustBeNumericOrLogical} = []
    opts.dotSize (1,1){mustBeNumeric} = 55;
    opts.edmund_lines (1,1) {mustBeNumericOrLogical} = false;
    opts.edmund_flip (1,1) {mustBeNumericOrLogical} = false;
    opts.chi_thresh  (1,1){mustBeNumeric} = 0.0;
    opts.tburn  (1,1){mustBeNumeric} = 0.0;
    opts.checker (1,1) {mustBeNumericOrLogical}  = false;
    opts.period (:,:) {mustBeNumericOrLogical} = true;
    opts.peaks (:,:)  = [];
    opts.dim (:,:) = [  600         400]
    opts.pad (1,1) {mustBeNumericOrLogical}  =0.05;
    opts.cycle_thresh (1,1) {mustBeNumeric} = 1.5;
    opts.raster (1,1) {mustBeNumericOrLogical} = false;
    opts.title (1,1)  = true;
    opts.xlabel (1,:) {mustBeTextScalar} = 'Oscillator #';
    opts.labelsize (1,1) {mustBeNumeric} = 18;
end
show_period = length(opts.period)>1 || false;
if ~isempty(opts.dim)
    opts.dim(2) =  opts.dim(2)+100*any(any(opts.period~=false));
end

upper_style = {'XTick',[], 'FontSize',14};
lower_style = { 'FontSize',opts.labelsize};

if ~isempty(sol)

    N=size(sol,2);
    if N>1
        mid=ceil(N/2)+(1-mod(N,2))/2;
        oscs = (1:N)-mid;
    else
        mid=1;
        oscs=[0];
    end

    peaks_mid = cell2mat(get_sol_peaks(sol(:,floor(mid)), tmesh, MinPeakDistance=opts.MinPeakDistance));
    if length(opts.period)>1
        T_mid = mean(diff(peaks_mid(peaks_mid >= opts.period(1,1) & peaks_mid< opts.period(1,2))));
    else
        T_mid = mean(diff(peaks_mid(peaks_mid >= opts.tmin & peaks_mid<opts.tmax)));
    end


    % extra.tmin = tmesh(1);



    tmin = max(opts.tmin, tmesh(1));
    tmax = min(tmesh(end), opts.tmax);


else

    peaks=opts.peaks;
    tmesh =[peaks{:}]; %this is just a temp variable
    tmesh = min(tmesh):max(tmesh);
    N=length(peaks);
    mid=ceil(N/2)+(1-mod(N,2))/2;
    oscs = (1:N)-mid;
    tmin =  max(opts.tmin, min(cellfun(@(x) min(x),peaks)));
    peaks_mid = peaks{mid};
    tmax=min(max(cellfun(@(x) max(x),peaks)),opts.tmax);
    T_mid = mean(diff(peaks_mid(peaks_mid >= opts.tmin & peaks_mid<opts.tmax)));

end




if ~isempty(sol)
    imin = find(tmesh >= tmin-4*T_mid,1);
    imax = find(tmesh <= tmax+4*T_mid,1, "last");
    if length(opts.period)>1
        sol0=sol;
        tmesh0=tmesh;
    end
    sol=sol(imin:imax,:);
    tmesh=tmesh(imin:imax);
    peaks = get_sol_peaks(sol, tmesh, MinPeakDistance=opts.MinPeakDistance, MinPeakProminence = opts.MinPeakProminence, MinPeakHeight=opts.MinPeakHeight);
else
    peaks = cellfun(@(x) x(x>= tmin-4*T_mid & x<= tmax+4*T_mid ), peaks,'UniformOutput', false);
end


if length(opts.period)>1

    peak_times_slice = cell(size(opts.period,1),1);
    mean_period_slice = cell(size(opts.period,1),1);
    period_err_slice = cell(size(opts.period,1),1);

    for i=1:size(opts.period,1)
        inds = tmesh0>=opts.period(i,1)&tmesh0<=opts.period(i,2);
        peaks_slice = get_sol_peaks(sol0(inds,:), tmesh0(inds), MinPeakDistance=opts.MinPeakDistance, MinPeakProminence = opts.MinPeakProminence, MinPeakHeight=opts.MinPeakHeight);
        peak_times_slice{i} = cellfun(@(p) p, peaks_slice,'UniformOutput', false);
        mean_period_slice{i} = cellfun(@(t) mean(diff(t)), peak_times_slice{i});
        period_err_slice{i} = cellfun(@(t) std(diff(t)), peak_times_slice{i});
    end

end

peak_times = cellfun(@(p) p, peaks,'UniformOutput', false);
mean_period = cellfun(@(t) mean(diff(t)), peak_times);
period_err = cellfun(@(t) std(diff(t)), peak_times);


if length(opts.period)==1


    % T_mid = mean(diff(tmid));
    T_mid = mean_period(floor(mid));
end
mean_period = mean_period/T_mid;
period_err = period_err/T_mid;

[ovec,tvec , phi_l, phi_r, phi_ls, phi_rs] = get_phases(peaks);
[chi, chis] = phase_asymmetry(phi_l, phi_r, phi_ls, phi_rs);
tmid = tvec(ovec==floor(mid)); % peaks_mid(peaks_mid>tmesh(1)&peaks_mid<tmesh(end));%

cycle0 = ceil((tmid(1)-tmin)/T_mid);


tiledlayout(2+any(any(opts.period~=false)), 1);

if any(any(opts.period~=false))
    nexttile(1);
    colororder(["k","r","c"])
    if length(opts.period)>1
        hold on
        labels = cell(size(opts.period,1),1);
        for i=1:size(opts.period,1)
            errorbar(oscs, mean_period_slice{i}/T_mid, period_err_slice{i}/T_mid);
            ci=round(min(cellfun(@ min,peak_times_slice{i}))/T_mid)+cycle0;
            cf = round(max(cellfun(@ max,peak_times_slice{i}))/T_mid)+cycle0;
            labels{i}=['Cycles ' num2str(ci) '-' num2str(cf)];
        end
        mean_period = cell2mat(mean_period_slice)/T_mid;
        period_err = cell2mat(period_err_slice)/T_mid;

        hold off
        hleg=legend(labels,'Location','Best');
        set(hleg,'fontsize',10)
    else
        errorbar(oscs, mean_period, period_err);
    end

    % set(gca, upper_style{:})

    range_err = max(mean_period+period_err)-min(mean_period-period_err);
    range = max(mean_period)-min(mean_period);
    % mean_err = mean(period_err);
    if 2*range<range_err
        if range_err>1e-8
            ylim([min(mean_period)-range_err*10, max(mean_period)+range_err*10])
        else
            set(gca,'YTick',[1])
        end

    else
        % ylim([min(mean_period)-range_err*0.1, max(mean_period)+range_err*0.1])
    end


    title('Period (normalized)','FontSize', 16)
    set(gca,  lower_style{:})
end
% colormap(gca, 'parula');
% ylabel('Period (norm.)', 'FontSize',14)
% xlabel('Oscillator #');





if opts.checker
    [psi, ~] = phase_checker(phi_l, phi_r, phi_ls, phi_rs);
    heatmap_vals = psi;
else
    heatmap_vals = chi;
end

xlim([0.5,N+0.5]-mid)


ax2 = nexttile(1+any(any(opts.period~=false)), [2,1]);
if isempty(tmid)
    return
end
if isempty(opts.cycle_min)
    cycle_min = floor((max(tmin, min(tvec))+opts.tburn)/T_mid);
    cycle_min=0;
else
    cycle_min=opts.cycle_min;
end



[cycle_map,~] = peaktimes_to_cycles(peaks, tmesh, tmin=opts.tmin, tburn=opts.tburn, tmax=opts.tmax , ind=floor(mid));


cycle = cycle_map(tvec);
if isempty(opts.cycle_max)
    cycle_max = cycle_map(min(tmax, max(tvec)));
else
    cycle_max=opts.cycle_max;
end
pad=opts.pad;
range=(max(cycle)-min(cycle));

yrange = [cycle_min-pad*range cycle_max+(pad)*range];

if opts.border && ~opts.raster
    % hold on;
    scatter(ovec-mid, cycle, opts.dotSize, [1, 1, 1]*0.6);
    hold on;
end

if ~opts.raster
    scatter(ovec-mid, cycle, opts.dotSize, heatmap_vals, "filled");
else
    inds = ~isnan(heatmap_vals);
    [xg,yg] = meshgrid(min(ovec-mid):max(ovec-mid), cycle_min:cycle_max);
    interp = scatteredInterpolant(ovec(inds)'-mid, cycle(inds)', heatmap_vals(inds)');
    s=pcolor(xg,yg,interp(xg,yg));
    set(s,'EdgeColor','none');
end

if opts.border
    hold off;
end


if ~isempty(opts.xlines)
    hold on;
    for i=opts.xlines
        h=plot([i, i], yrange ,'-k');
        uistack(h,"bottom");
    end
    hold off;

end

xlabel(opts.xlabel,'FontSize', opts.labelsize);
ylabel('Cycle #','FontSize',opts.labelsize)

set(gca, lower_style{:})


xlim([0.5,N+0.5]-mid)


ylim(yrange)
if ~opts.checker
    if islogical(opts.title) && opts.title==true
        opts.title='Phase Asymmetry';
    end
    map = customcolormap_preset('purple-white-green');
    map=sqrt(map);
    colormap(gca, map);
    h=colorbar();
    if ~isempty(opts.clim)
        clim(opts.clim)
    end
    title(h, '\chi','FontSize',opts.labelsize);
else
    if islogical(opts.title) && opts.title==true
        opts.title='Anti-Phase Metric';
    end
    colormap(gca, hot);
    clim([0,1])
    h=colorbar();
    title(h, '\psi','FontSize',opts.labelsize);
end

if ischar(opts.title)
    title(ax2, opts.title,'FontSize', 16)
end

if ~isempty(opts.dim)
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [ pos(1:2)        opts.dim] )
end

if ~isempty(opts.phi_scatter)
    figure(opts.phi_scatter);
    phi_scatter(cell2mat([phi_rs{:}]), cell2mat([phi_ls{:}]));
end

%Edmund style wave lines
if opts.edmund_lines && ~ opts.raster
    if opts.edmund_flip
        flip=-1;
    else
        flip=1;
    end
    x = {};
    y = {};
    for n=1:N
        for j=1:length(chis{n})
            if flip*chis{n}(j)>flip*opts.chi_thresh && n>=3
                next = find(peaks{n-2}>peaks{n}(j),1);
                if ~isempty(next) && flip*chis{n-2}(next)>=flip*opts.chi_thresh
                    y0=cycle_map(peaks{n}(j));
                    yf=cycle_map(peaks{n-2}(next)) ;
                    if yf-y0<opts.cycle_thresh
                        x{end+1} = [n,n-2]-mid;
                        y{end+1} = [y0,yf];
                    end
                    % line([n,n-2]-mid,[y0,yf],'color','k');
                end
            elseif flip*chis{n}(j)<-flip*opts.chi_thresh && n<=N-2
                next = find(peaks{n+2}>peaks{n}(j),1);
                if ~isempty(next) && flip*chis{n+2}(next)<=-flip*opts.chi_thresh
                    y0= cycle_map(peaks{n}(j));
                    yf= cycle_map(peaks{n+2}(next)) ;
                    if yf-y0<opts.cycle_thresh
                        x{end+1} = [n,n+2]-mid;
                        y{end+1} = [y0,yf];
                    end
                    % line([n,n+2]-mid,[y0,yf],'color','k');
                end
            end
        end
    end
    x=cell2mat(x');
    y=cell2mat(y');
    h=line(x',y','color',[1 1 1]*0.4);
    uistack(h,"bottom");

    % x = {};
    % y = {};
    % for n=1:N
    %     for j=1:length(chis{n})
    %         if chis{n}(j)<-extra.chi_thresh && n>=3 && n~=mid+1
    %             next = find(peaks{n-2}>peaks{n}(j),1);
    %             if ~isempty(next) && chis{n-2}(next)<-extra.chi_thresh
    %                 y0=cycle_min + interp1(tmid,cycle_inds,peaks{n}(j),'linear','extrap');
    %                 yf=cycle_min + interp1(tmid,cycle_inds,peaks{n-2}(next),'linear','extrap');
    %                 x{end+1} = [n,n-2]-mid;
    %                 y{end+1} = [y0,yf];
    %                 % line([n,n-2]-mid,[y0,yf],'color','k');
    %             end
    %         elseif chis{n}(j)>extra.chi_thresh && n<=N-2  && n~=mid-1
    %             next = find(peaks{n+2}>peaks{n}(j),1);
    %             if ~isempty(next) && chis{n+2}(next)>extra.chi_thresh
    %                 y0=cycle_min + interp1(tmid,cycle_inds,peaks{n}(j),'linear','extrap');
    %                 yf=cycle_min + interp1(tmid,cycle_inds,peaks{n+2}(next),'linear','extrap');
    %
    %                 x{end+1} = [n,n+2]-mid;
    %                 y{end+1} = [y0,yf];
    %                 % line([n,n+2]-mid,[y0,yf],'color','k');
    %             end
    %         end
    %     end
    % end
    % x=cell2mat(x');
    % y=cell2mat(y');
    % h=line(x',y','color',[1 1 1]*0);
    % uistack(h,"bottom");
end


end