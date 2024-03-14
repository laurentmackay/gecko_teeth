figure(1)
tspan = [0,15e4];  % start and end times
n=5;
% v0 =  0; w0 = 0; % initial values
IC = [repmat( 0 ,n,1) repmat( 0 , n,1)];
IC=randn(n,2);
% IC=[ -0.4110   0.6636 -.4110    -0.2258 0.2128 -.22]
%N=5 phase asymmetry
% IC=[0.0496   -0.2220    0.0003   -0.2220    0.0496   -0.0043    0.0976   -0.0002    0.0976   -0.0043];

% IC(3)=randn()*3;
% IC(ceil(n/2))=1.8;
% IC(ceil(n/2),2)=-.1
% Call ode45

%%
tspan = [0,6e5];

options = odeset('MaxStep',5, BDF='on');
D=0.1;
b=0.105;
fun  = @ SCH_zf_array;
IC(:)
[t, vw] = ode15s( fun(n , D , b) , tspan, IC(:), options);
% Extract individual solution values
v = vw(:,1:n);
w = vw(:,n+1:end);
% % Plot results
% plot(t,v,'r',t,w,'b'),grid
% xlabel('t'),ylabel('v and w')
% legend('v','w')

figure(1)
t_disp=1e4;
i_disp=find(t>tspan(end)-t_disp,1);
h=imagesc(1:n,t(i_disp:end)-t(i_disp),v(i_disp:end,:));
% h=pcolor(1:n,t,v);
% set(h, 'EdgeColor', 'none');
colorbar
IC = vw(end,:);
disp('done')

%%
figure(2)
per = zeros(n,1);
for i=1:n
    vmax =  max(v(i_disp:end,i));
    vmin =  min(v(i_disp:end,i));
    range = vmax-vmin;
    findpeaks(v(i_disp:end,i),'MinPeakProminence',0.5*range)
[~,p] = findpeaks(v(i_disp:end,i),'MinPeakProminence',0.5*range);
per(i) = mean(diff(t(i_disp+(p-1))));
drawnow
end
per

%%
NTST=1000
T=mean(per)

tspan = linspace(0,T, NTST+1);
% options = odeset('MaxStep',T);
[tout, vwout] = ode45( fun(n , D, b) ,[0,T], vw(end,:));

rel_delta = (vwout(1,:)- vwout(end,:))./ vwout(1,:)
vwout(end,:)=vwout(1,:);
return
%%
fid = fopen('/home/lmackay/atest/five_boy_dead.dat','w','n')
fprintf(fid, [ repmat('%2.9E ', 1, 1+2*n) '\n'], [tout(1:end) vwout(1:end,:)]')


