function plot_switching_times(N0, T0, theta0, gamma, opts)
arguments
N0 (1,1) {mustBeNumeric}
T0 (1,1) {mustBePositive}
theta0 (1,:) {mustBeNumeric}
gamma (1,1) {mustBeNumeric}
opts.labelsize (1,1) {mustBePositive} = 14
opts.bottom (1,1) {mustBeNumericOrLogical} = false
opts.k (1,:) {mustBeNumeric} = 0:10
opts.tburn (1,1) {mustBeNumeric} = 0
end
hold on
N=2*N0+1;
n=-N0:N0;

dr=mod([theta0(1:end-1)-theta0(2:end), nan],2*pi)/(2*pi);
dl=mod([nan,theta0(2:end)-theta0(1:end-1)],2*pi)/(2*pi);

yrange=diff(ylim);
pos=get(gca,'Position');

c0=opts.tburn/T0;
% c0 = floor(c0)- mod(c0,1);
for k=opts.k


   
    for kk=(2*k):(2*k+1)
  
    swc=[(dl-dr+kk*sign(n)).*(1./(N0*gamma*2*sign(n))).*(N0+gamma*abs(n-1)).*(N0+gamma*abs(1+n))] + c0;
    h=plot(n, swc,':', 'LineWidth', 2, 'Color',[0,0,0,0.9]);
    angle=180*atan2(diff(pos(4)*swc(end-3:end-2))/yrange,pos(3)/N)/pi;
    text(n(end)-0.5, swc(end-1)+(0.02)*yrange, ['k=' num2str(kk)],'HorizontalAlignment','right','VerticalAlignment','bottom','BackgroundColor',[1,1,1, 0.85],'Clipping','on','Margin',0.5,'Rotation',angle,'FontSize',opts.labelsize)
    if opts.bottom
        uistack(h,"bottom")
    end
    % u
    end

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
