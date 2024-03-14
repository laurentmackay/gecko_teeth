function q=PlotPeriodPhasePhaseAsym(EruptionTimes,dotSize,contextShort,contextLong,saveEPS)

% Font sizes for plots
FS1 = 28;  % Panel letters A, B, C+b(1,
FS2 = 20;  % Panel titles and axes label
FS3 = 16;  % Color bar label

EndDay = max(EruptionTimes(:,2));
LetterLocation = EndDay*1.2;

JawEndL = -37.5; %min(EruptionTimes(:,1)) - 2;
JawEndR = 37.5; %max(EruptionTimes(:,1)) + 2;

%% Plot period
    lapsedTime=[-1; EruptionTimes(2:end,2)-EruptionTimes(1:end-1,2)];
    indDefined=find(lapsedTime>=0);
    indUndefined=find(lapsedTime<0);
    count=0;
    for k=min(EruptionTimes(indDefined,1)):max(EruptionTimes(indDefined,1))
        count=count+1;
        indToothk=find(EruptionTimes(indDefined,1)==k);
        MeanPeriod(count)=mean(lapsedTime(indDefined(indToothk)));
        StDevPeriod(count)=std(lapsedTime(indDefined(indToothk)));
        
    end
        [EruptionTimes(indDefined,1) lapsedTime(indDefined)]
        [MeanPeriod' [min(EruptionTimes(indDefined,1)):max(EruptionTimes(indDefined,1))]']

    mn=mean(lapsedTime(indDefined));
    st=std(lapsedTime(indDefined));
    
    f=figure(1);
    f.Position = [1101 180 780 841];
%    f.Position = [  1017 224 616 761]; 

    ax(1)=subplot(3,1,1);
    cla
    hold on
    text(-50,LetterLocation,'A','fontsize',FS1)

    if length(contextShort)==5
        if contextShort=='LG190'
            
            coordsInt =    [50     7;    66   249;    68   249;    53    14;
                            70   245;    55    11;    57     4;    72   245;
                            75   256;    60    11;    77   249;    62     4;
                            26     7;    10   256;    13   245;    28    11;
                            30    14;    15   252;    17   259;    33     4;
                            35     7;    20   249;    37    14;    22   256;
                            ];
%         coordsInt =    [
%                         11     4;     8    46;     8    84;    13    11;    16     4;     8   126;
%                          8   165;    18    11;    21     4;     8   207;     8   252;    23    11;
%                         26     7;    10   256;    13   245;    28    11;    30    14;    15   252;
%                         17   259;    33     4;    35     7;    20   249;
%                         37    14;    22   256;    25   245;    40     7;
%                         42    11;    27   249;    29   256;    44    18;    44    49;    31   259;
%                         34   245;    44    84;    44   123;    36   256;
%                         39   249;    44   161;    44   200;    41   256;
%                         44   242;    45   259;
%                         44   200;    47   256;    44   161;    49   252;
%                         44   123;    51   245;    54   256;    44    84;    44    49;
%                         56   249;    59   256;    44    18;    
%                         46    14;    61   249;
%                         64   256;    48     7;    
%                         50     7;    66   249;    68   249;
%                         53    14;    70   245;    55    11;    57     4;    72   245;
%                         75   256;    60    11;    77   249;    62     4;
%                         79   242;
%                         65    21;    67    11;    80   214;    81   193;    69     4;
%                         72    11;    81   147;    80    91;    74     4;    77    11;
%                         80    49
%                         ];

            coordsInt(:,1) = coordsInt(:,1) - 44;
            for k=1:2:size(coordsInt,1)
                plot([coordsInt(k,1) coordsInt(k+1,1)],[coordsInt(k,2) coordsInt(k+1,2)],'color',[0.6 0.6 0.6]) %[0 0.4470 0.7410]
            end
            
            coordsInt = [
                        0        84;
                        18       35;
                        0        84;
                       -18       39;
                        0       123;
                      -28        39;
                        0       123;
                       30        39;
                       -1       102;
                      -23        39;
                        1       102;
                       25        39];
            for k=1:2:size(coordsInt,1)
                plot([coordsInt(k,1) coordsInt(k+1,1)],[coordsInt(k,2) coordsInt(k+1,2)],'--','color',[0.6 0.6 0.6]) %[0 0.4470 0.7410]
            end
            

        end
    end
    
    
        if length(contextShort)==5
            if contextShort=='LG173' | contextShort=='LG174' | contextShort=='LG178' | contextShort=='LG219'
                EruptionTimes(:,1) = -EruptionTimes(:,1);
            end
        end
        
    scatter(EruptionTimes(indDefined,1),EruptionTimes(indDefined,2),dotSize,lapsedTime(indDefined),'filled')
    % scatter(EruptionTimes(indUndefined,1),EruptionTimes(indUndefined,2),dotSize,[0.6 0.6 0.6])
    % xlim([0 81]);
%    xlabel('Tooth Position','interpreter','latex','FontSize',42)
%    ylabel('Day of eruption','interpreter','latex','FontSize',FS2)
    ylabel('Day of eruption','FontSize',FS2)
    mycolormap1=jet;
    mycolormap1=mycolormap1(1:128,:);
    colormap(ax(1),mycolormap1);

    colorbar;
    h=colorbar;
    caxis([0 65])
    if length(contextShort)==8
        if contextShort=='CrocData'
            caxis([250 650])
        end
    end
    ylabel(h, 'Days since last eruption','FontSize',FS3)
    title({contextLong,'Days since last eruption'},'FontSize',FS2)
    ax = gca;
    ax.FontSize = FS2; 
%    ax.TickLabelInterpreter = 'latex';
    h.FontSize = FS3; 
%    h.TickLabelInterpreter = 'latex';
    
    axis([JawEndL JawEndR -15 EndDay+15])
%        axis([JawEndL JawEndR -15 232])

%     if saveEPS
%         filename=sprintf('Figures/%s/period_%s.eps',contextShort,contextShort)
%         saveas(gcf,filename,'epsc')
%     end
    
%% Best fit to T0+a*N+b*t

N = EruptionTimes(indDefined,1);
T = EruptionTimes(indDefined,2);
P = lapsedTime(indDefined);
% indKeep = find(N>=0);
indRightSide = (N>=0);
indSurgicalSite = (N>14).*(N<23).*(T>49).*(T<162);
indKeep = find(indRightSide.*(1-indSurgicalSite));
N = N(indKeep);
T = T(indKeep);
P = P(indKeep);
X = [ones(size(N)) N T];
%[b,bint] = regress(P,X);

% Fit for left side of LG190
% b =
% 
%    30.1144
%    -0.0761
%     0.0455
% 
% 
% bint =
% 
%    29.2731   30.9557
%    -0.1035   -0.0487
%     0.0411    0.0500

% For the right side of LG190
% b =
% 
%    33.4673
%     0.1188
%     0.0294
% 
% 
% bint =
% 
%    29.4126   37.5221
%    -0.0096    0.2471
%     0.0081    0.0508

% For the right side of LG190 with just-post surgical eruption omitted
% b =
% 
%    30.6069
%     0.1059
%     0.0392
% 
% 
% bint =
% 
%    28.6511   32.5626
%     0.0446    0.1671
%     0.0289    0.0494

% figure(5)
% plot3(N,T,P)

%% Plot phase

    for k=1:size(EruptionTimes,1)
        % phase of current tooth relative to the previous and next tooth in the position to the left 
        indLeft = find(EruptionTimes(:,1) == EruptionTimes(k,1)-1); % returns indices of all teeth in the position to the left of tooth k
        indPrev = find(EruptionTimes(indLeft,2)<EruptionTimes(k,2));
        prev = max(EruptionTimes(indLeft(indPrev),2));
        indNext = find(EruptionTimes(indLeft,2)>=EruptionTimes(k,2));
        next = min(EruptionTimes(indLeft(indNext),2));
        leftPhase = (EruptionTimes(k,2)-prev) / (next-prev);
        % phase of current tooth relative to the previous and next tooth in the position to the right 
        indRight = find(EruptionTimes(:,1) == EruptionTimes(k,1)+1); % returns indices of all teeth in the position to the left of tooth k
        indPrev = find(EruptionTimes(indRight,2)<EruptionTimes(k,2));
        prev = max(EruptionTimes(indRight(indPrev),2));
        indNext = find(EruptionTimes(indRight,2)>=EruptionTimes(k,2));
        next = min(EruptionTimes(indRight(indNext),2));
        rightPhase = (EruptionTimes(k,2)-prev) / (next-prev);
        if isempty(leftPhase)|isempty(rightPhase)
            avgPhase(k) = NaN;
            asymPhase(k) = NaN;
        else
            avgPhase(k) = (rightPhase+leftPhase)/2;
            asymPhase(k) = rightPhase-leftPhase;
        end
    end
    
    ind=find(asymPhase>0);
    mean(asymPhase(ind));
    ind=find(asymPhase<0);
    mean(asymPhase(ind));

    indNaN=find(isnan(avgPhase));
    indNotNaN=find(~isnan(avgPhase));

    test=hot(256);
    test=[ flipud(test(50:177,:));test(50:177,:)];

    ax(2)=subplot(3,1,2);
    cla
    hold on
    
    text(-50,LetterLocation,'B','fontsize',FS1)

    scatter(EruptionTimes(indNotNaN,1),EruptionTimes(indNotNaN,2),dotSize,avgPhase(indNotNaN),'filled')
    % scatter(EruptionTimes(indNaN,1),EruptionTimes(indNaN,2),dotSize,[0.6 0.6 0.6])
    %scatter(xvalues(phase_vector_r(:,2)>1.3),yvalues(phase_vector_r(:,2)>1.3),150,'k')
    % xlim([0 81]);
%    xlabel('Tooth Position','interpreter','latex','FontSize',42)
    ylabel('Day of eruption','FontSize',FS2)
    colormap(ax(2),test)
    h=colorbar;
%    ylabel(h, 'Average relative phase','FontSize',FS2)
    ylabel(h, '$\displaystyle\frac{\phi_r+\phi_l}{2}$','interpreter','latex','FontSize',FS3*1.3)
    % caxis([0 2*pi])
    % caxis([0 0.75])
    caxis([0 1])
    title(['Average relative phase'],'FontSize',FS2)
    ax = gca;
    ax.FontSize = FS2; 
%    ax.TickLabelInterpreter = 'latex';
%    h.FontSize = FS3; 
%    h.TickLabelInterpreter = 'latex';
    h.Ticks = linspace(0, 1, 5) ; %Create 8 ticks from zero to 1
    h.TickLabels = num2cell(0:0.25:1) ; 
    axis([JawEndL JawEndR -15 EndDay+15])
%    axis([JawEndL JawEndR -15 232])
%     if saveEPS
%         filename=sprintf('Figures/%s/avg_phase_%s.eps',contextShort,contextShort);
%         saveas(gcf,filename,'epsc')
%     end
    
%% Plot phase asymmetry    

    ax(3)=subplot(3,1,3);
    cla
    hold on
    
    text(-50,LetterLocation,'C','fontsize',FS1)

    scatter(EruptionTimes(indNotNaN,1),EruptionTimes(indNotNaN,2),dotSize,asymPhase(indNotNaN),'filled')
    % scatter(EruptionTimes(indNaN,1),EruptionTimes(indNaN,2),dotSize,[0.6 0.6 0.6])
    % xlim([0 81]);
    xlabel('Tooth Position','FontSize',FS2)
    ylabel('Day of eruption','FontSize',FS2)
    colormap(ax(3),sqrt(customcolormap_preset('purple-white-green')));
    h=colorbar;
%    ylabel(h, 'Asymmetry of relative phase','FontSize',FS3)
    ylabel(h, '$\displaystyle \phi_r-\phi_l$','interpreter','latex','FontSize',FS3*1.3)
    % caxis([0 2*pi])
    caxis(2/3*[-0.3 0.3])
    title(['Asymmetry of relative phase'],'FontSize',FS2)
    ax = gca;
    ax.FontSize = FS2; 
%    ax.TickLabelInterpreter = 'latex';
%    h.FontSize = FS3; 
%    h.TickLabelInterpreter = 'latex';
    h.Ticks = linspace(-0.3, 0.3, 5) ; %Create 8 ticks from zero to 1
    h.TickLabels = num2cell(-0.3:0.15:0.3) ; 

    axis([JawEndL JawEndR -15 EndDay+15])
%    axis([JawEndL JawEndR -15 232])
%     figure(3)
%     clf
%     hold on
%     scatter(EruptionTimes(indNotNaN,1),EruptionTimes(indNotNaN,2),dotSize,asymPhase(indNotNaN),'filled','MarkerEdgeColor',[0.8 0.8 0.8]);hold on
%     scatter(EruptionTimes(indNaN,1),EruptionTimes(indNaN,2),dotSize,[0.6 0.6 0.6]);hold on
%     % xlim([0 81]);
%     xlabel('Tooth Position','interpreter','latex','FontSize',42)
%     ylabel('Day of eruption','interpreter','latex','FontSize',42)
%     colormap(gca,sqrt(customcolormap_preset('purple-white-green')));
%     h=colorbar;
%     ylabel(h, 'Asymmetry of relative phase','interpreter','latex','FontSize',42)
%     % caxis([0 2*pi])
%     caxis(2/3*[-0.3 0.3])
%     title(['Asymmetry of relative phase'],'interpreter','latex','FontSize',42)
%     ax = gca;
%     ax.FontSize = 30; 
%     ax.TickLabelInterpreter = 'latex';
%     h.FontSize = 30; 
%     h.TickLabelInterpreter = 'latex';
%     h.Ticks = linspace(-0.3, 0.3, 5) ; %Create 8 ticks from zero to 1
%     h.TickLabels = num2cell(-0.3:0.15:0.3) ; 

    
%% Save plot to eps
    if saveEPS
        filename=sprintf(['Figures/Fig_%s.eps'],contextShort);
        saveas(gcf,filename,'epsc')
    end
    
    
end