function [ls, gofs] = fit_windowed_periods(peaks, tmesh, wins, sym, opts)
arguments 
peaks (:,:) {mustBeA(peaks, 'cell')}
tmesh (1,:) {mustBeNumeric}
wins (:,2) {mustBeNumeric}
sym (1,1) {mustBeNumericOrLogical} = true
opts.mask (1,:) {mustBeNumericOrLogical} = true(1,length(peaks));
opts.m (1,1) {mustBeNumeric} = 1;
opts.plot (1,1) {mustBeNumericOrLogical} =false;
end

inds=1:length(peaks);
mid=ceil(length(peaks)/2);

periods = cellfun(@(t) diff(tmesh(t)), peaks,'UniformOutput', false);

mask=opts.mask;
m=opts.m;


mean_period_win = arrayfun(@(start, stop) cellfun(@(p,t) mean(p(tmesh(t(2:end))>=start & tmesh(t(2:end))<=stop)),periods, peaks),wins(:,1), wins(:,2) , 'UniformOutput', false);
mean_period_win = cell2mat(mean_period_win);
mean_period_win(:,~mask)=nan;
% subplot(1,2,1)


ls={};gofs={};


func = @(p0, gamma_l, gamma_r, x) p0*(gamma_l*(abs((x-mid)/(mid-inds(1))).^m).*(x<=mid)+...
    gamma_r*(abs((x-mid)/(inds(end)-mid)).^m).*(x>=mid)+1);
ft = fittype(func, coefficients=["p0" "gamma_l" "gamma_r"],independent="x");

func_sym = @(T0, DeltaT, x) T0+DeltaT*((abs((x-mid)/(mid-inds(1))).^m).*(x<=mid)+(abs((x-mid)/(inds(end)-mid)).^m).*(x>=mid));
ft_sym = fittype(func_sym, coefficients=["T0" "DeltaT"],independent="x");


for i = 1:size(wins,1)
    mask_2 = mask & ~isnan(mean_period_win(i,:));
    % good_2 = ~isnan(mean_period_win(i,:));
    T0 = mean_period_win(i,inds==mid);
    TR = max(mean_period_win(i,inds>mid));
    TL = max(mean_period_win(i,inds<mid));
    if sym
        [ls{end+1}, gofs{end+1} ]= fit(inds(mask_2)',  mean_period_win(i,mask_2)',ft_sym, startpoint=[T0,(TR+TL-2*T0)/2]);
    else
        [ls{end+1}, gofs{end+1} ] = fit(inds(mask_2)',  mean_period_win(i,mask_2)',ft, startpoint=[mean_period_win(1,inds==mid),0.1,0.1]);
    end
end

if opts.plot
    times=mean(wins,2);
    figure(opts.plot);
    surf(times, inds-mid, mean_period_win');
    
    hold on
    surf(times, inds-mid, cell2mat(cellfun(@(cfit) cfit(inds),ls,'UniformOutput',false)),'EdgeColor','none','FaceColor',[0,0,0])
    title(['mean r^2=' num2str(mean(cellfun(@(gof) gof.rsquare, gofs))) ])
    hold off
end

end