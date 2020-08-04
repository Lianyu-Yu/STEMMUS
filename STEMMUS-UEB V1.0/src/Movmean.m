%% Energy balance CHK

%% bar figure
CK11=[51.8200000000000,53.6500000000000;75.3800000000000,69.4200000000000;4.35000000000000,3.86000000000000;2.57000000000000,2.47000000000000;-34.6900000000000,-29.9500000000000;-0.550000000000000,-0.550000000000000]
% CK11=[53.2700000000000,54.7500000000000,57.3600000000000;72.1800000000000,61.0100000000000,62.3300000000000;4.15000000000000,3.86000000000000,4.99000000000000;2.47000000000000,2.50000000000000,2.39000000000000;-32.3900000000000,-22.4600000000000,-27.6000000000000;-0.320000000000000,-0.320000000000000,-0.340000000000000]
figure
b=bar(CK11(1:6,:)',0.5,'stacked')
% hold on
% bar(CK11(5,:)',0.5,'g')
% bar(CK11(6,:)',0.5,'r')
ck=3
for i=1:ck
%     if k==1
        y(1,i)=CK11(1,i)/2;
        text(i,y(1,i),num2str(CK11(1,i)),'HorizontalAlignment','center','VerticalAlignment','middle')
%     end
    for k=2:6
        y(k,i)=CK11(k,i)/2+sum(CK11(1:k-1,i));
%         if k==5
%          y(k,i)=y(k,i)+5;
%         end
        text(i,y(k,i),num2str(CK11(k,i)),'HorizontalAlignment','center','VerticalAlignment','middle')
    end
%     if k==4
%         y(k,i)=CK11(k,i)/2+sum(CK11(1:k-1,i))+10;
%         text(i,y(k,i),num2str(CK11(k,i)),'HorizontalAlignment','center','VerticalAlignment','middle')
%     end
end
% for i=1:ck
% %     if k==1
%         y(5,i)=CK11(5,i)/2;
%         text(i,y(5,i),num2str(CK11(5,i)),'HorizontalAlignment','center','VerticalAlignment','middle')
% %     end
% %     for k=2:4
% %         y(k,i)=CK11(k,i)/2+sum(CK11(1:k-1,i));
% %         text(i,y(k,i),num2str(CK11(k,i)),'HorizontalAlignment','center','VerticalAlignment','middle')
% %         
% %     end
% end

% for i=1:ck
% %     if k==1
%         y(6,i)=CK11(6,i)/2;%+CK11(5,i);
%         text(i,y(6,i),num2str(CK11(6,i)),'HorizontalAlignment','center','VerticalAlignment','middle')
% %     end
% %     for k=2:4
% %         y(k,i)=CK11(k,i)/2+sum(CK11(1:k-1,i));
% %         text(i,y(k,i),num2str(CK11(k,i)),'HorizontalAlignment','center','VerticalAlignment','middle')
% %         
% %     end
% end
str = {'unCPLD';'unCPLD-FT'; 'CPLD'};%; 'CPLD-NWB'
% bar(x)
set(gca,'ylim',[0,100]);  %[0 10 20 30 40 50 60 70 80 90]
% set(gca,'xlim',[0,1]);  %[0 10 20 30 40 50 60 70 80 90]
% set(gca,'xtick',0:0.2:1,'xticklabel',{0 0.2 0.4 0.6 0.8 1 })
set(gca,'ytick',0:20:100,'yticklabel',{0 20 40 60 80 100 })

set(gca, 'XTickLabel',str, 'XTick',1:numel(str))
legend('Tv','Es','E_{IN}','E_{SN}','\Delta Vs','L_{K}','Location','northoutside','Orientation','horizontal')
legend boxoff
ylabel('Percentage (%)') % left y-axis 
hold on 
box off
hold off

figure(5)
%% fig 1
subplot(2,1,1)
x=1:1:length(Pr);
% QE12(QE12(:,5)==0)=nan;
ind=find(QE12(:,5)==0);
QE12(ind,5)=nan;
ind=find(QE12(:,6)==0);
QE12(ind,6)=nan;
[hAx,hLine1,hLine2]=plotyy(x,QE12(:,6),x,QE12(:,4));
hold on
y=QE12(:,5);
plot(y,'Color','k','linestyle','-','linewidth',1)
% hold on
hLine1.LineWidth = 1;
hLine2.LineWidth = 1;
% % set(gca,'xlim',[0 270])
set(gca,'xtick',1:90*24:length(Pr),'xticklabel',{'3/25' '6/23' '9/21' '12/20' '3/20' '6/18' '9/16' '12/15' '3/15' '6/13'})
% xlabel('Date (mm/dd)')
% % set(gca,'xtick',0:20:266,'xticklabel',{324 344 364 19 39 59 79 99 119 139 159 179 199 219})
% % xlabel('Day of Year(DoY)')
ylabel(hAx(1),'Cum. ET/P (mm)') % left y-axis 
ylabel(hAx(2),'P (mm)','position',[length(Pr)+900,10,0]) % right y-axis
% ylim(hAx(1),[0 1000])
% ylim(hAx(2),[0 40])
set(hAx(1),'xlim',[0 length(Pr)]);
set(hAx(2),'xlim',[0 length(Pr)]);

set(hAx(1),'ylim',[0 1000]);
set(hAx(2),'ylim',[0 100]);
legend('Cum. P','Cum. ET','P')
% title('(a)','Position', [10, 370, 0])
title('(a)','Position', [0, 1020, 0])
% set(gca,'fontsize',14);
set(hAx(1),'ytick',[0:200:1000]);
set(hAx(2),'ylim',[0 100]);
set(hAx(2),'ytick',[0:10:20]);
set(hAx(2),'Ydir','reverse');
legend boxoff 
% set(Ax(2),'Ydir','reverse');
% ax=hAx(2);
% ax.XAxisLocation='top';
% axis off

%% %%%%%%%%%%% freezing front drawing
% figure
subplot(2,1,2)

DEPTH1=[2.5 5 10 20 40 60 100];
% DEPTH=[0.25	0.5	1	2	3	4	5	6	7	8	9	10	12	14	16	18	20	22	24	26	28	30	32	34	36	38	40	45	50	55	60	70	80	100	120	140	160];
% contourf(0:3168,DEPTH,SIMTEMP',50,'Color','none');
% hold on
OBSTEMP=Tdp1(1:end,1:1:7);
% OBSTEMP(3500:4608,1:4)=-99;
OBSTEMP(OBSTEMP<-999)=nan;
% figure
[c1,h1]=contour(1:length(OBSTEMP(:,1)),DEPTH1,OBSTEMP',[0 0],'Color','k','linestyle','-','linewidth',2);
% xlim([0 5067/48])
% set(handles,'xtick',0:10:5067/48)
% axis([0,5067/48,-inf,inf])
set(gca,'Ydir','reverse');
ax=gca;
ax.XAxisLocation='bottom';
set(gca,'xtick',1:90*24:length(OBSTEMP(:,1)),'xticklabel',{'3/25' '6/23' '9/21' '12/20' '3/20' '6/18' '9/16' '12/15' '3/15' '6/13'})
xlabel('Date (mm/dd)')
ylabel('Soil depth (cm)')
title('(b)','Position', [0, -1, 0])
% for i=1:floor(length(QE)/24/90)
% plot([i*90*24,i*90*24],ylim,'k--'); % ??x=8,9,10,11???
% end

%% plot the coup/uncoup comparison fig9 of GPP
NofM=length(QE);
figure
load('GPP_UPDATE.mat')
ind=find(NPP_L11(:,3)<=0);
NPP_L11(ind,3)=NaN;
ind=find(NPP_L11(:,1)==0);
NPP_L11(ind,1)=NaN;
ind=find(NPP_L11(:,6)==0);
NPP_L11(ind,6)=NaN;
plot(NPP_L11(:,1),NPP_L11(:,2),'s-m','linewidth',1,'MarkerEdgeColor','m',...
'MarkerFaceColor','m')%
hold on
plot(NPP_L11(:,3),NPP_L11(:,6),'ok','linewidth',1,'MarkerEdgeColor','k',...
'MarkerFaceColor','k','Markersize',4)
plot(NPP_L11(:,3),NPP_L11(:,4),'-r','linewidth',1)
plot(NPP_L11(:,3),NPP_L11(:,5),'-g','linewidth',1)
plot(NPP_L11(:,3),NPP_L11(:,7),'-b','linewidth',1)


% plot(NPP_L12(:,4),NPP_L12(:,7),':r','linewidth',3)

set(gca,'ylim',[0,15]);  %[0 10 20 30 40 50 60 70 80 90]
set(gca,'xlim',[85,floor(NofM/24)+90]);  %[0 10 20 30 40 50 60 70 80 90]
clear title xlabel ylabel
set(gca,'xtick',85:90:(floor(NofM/24)+90),'xticklabel',{'3/25' '6/23' '9/21' '12/20' '3/20' '6/18' '9/16' '12/15' '3/15' '6/13'})
xlabel('Date (mm/dd)')
ylabel('GPP (g C m^{-2} d^{-1})')
legend('MODIS','OBS','unCPLD','unCPLD-FT','CPLD','Location','North','Orientation','horizontal','NumColumns',3)
legend boxoff
% text(410,7,'unCPLD: R^2 = 0.74, RMSE= 1.74 g C m^{-2} d^{-1}')
% text(410,6,'CPLD  : R^2 = 0.74, RMSE= 1.77 g C m^{-2} d^{-1}')
% text(0,8,'R^2 = 0.75')
% text(0,7,'RMSE= 0.71 g C m^{-2} d^{-1}')
title('(a)','Position', [85, 15.4, 0])

%% LAI SIMULATIONS FIG. 9
%%%%% FOR THE COUPLED SIMULATIONS
% ind=find(LAI_L1(:,1)==0);
% LAI_L1(ind,3)=NaN;
load('FIG8-LAI.mat')
figure
NofM=length(QE);
ind=find(LAI_L11(:,1)==0);
LAI_L11(ind,1)=NaN;
ind=find(LAI_L11(:,3)==0);
LAI_L11(ind,3)=NaN;
plot(LAI_L11(:,1),LAI_L11(:,2),'s-m','linewidth',1,'MarkerEdgeColor','m',...
'MarkerFaceColor','m')%
hold on
plot(LAI_L11(:,3),LAI_L11(:,4),'-R','linewidth',1)
plot(LAI_L11(:,3),LAI_L11(:,5),'-g','linewidth',1,'Markersize',8)
plot(LAI_L11(:,3),LAI_L11(:,6),'-b','linewidth',1,'Markersize',8)

set(gca,'ylim',[0,8]);  %[0 10 20 30 40 50 60 70 80 90]
set(gca,'xlim',[85,(floor(NofM/24)+90)]);  %[0 10 20 30 40 50 60 70 80 90]

set(gca,'xtick',85:90:(floor(NofM/24)+90),'xticklabel',{'3/25' '6/23' '9/21' '12/20' '3/20' '6/18' '9/16' '12/15' '3/15' '6/13'})
set(gca,'ytick',0:2:8,'yticklabel',{0 2 4 6 8})

xlabel('Date (mm/dd)')
ylabel('LAI (m^2/m^2)')
legend('MODIS','unCPLD','unCPLD-FT','CPLD','Location','North','Orientation','horizontal')
legend boxoff
title('(b)','Position', [85, 8.2, 0])

%% plot the coup/uncoup comparison fig of GPP
load('NEE_UPDATE.mat')
% ind=find(NEE1(:,1)<=0);
% NEE1(ind,1)=NaN;
ind=find(NEE1(:,1)==0);
NEE1(ind,1)=NaN;
ind=find(NEE1(:,4)==0);
NEE1(ind,4)=NaN;
% plot(NEE1(:,4),NEE1(:,1),'om','linewidth',1,'MarkerEdgeColor','m',...
% 'MarkerFaceColor','m')%
figure
plot(NEE1(:,4),NEE1(:,1),'ok','linewidth',1,'MarkerEdgeColor','k',...
'MarkerFaceColor','k','Markersize',4)
hold on
plot(NEE1(:,4),NEE1(:,2),'-r','linewidth',1)
plot(NEE1(:,4),NEE1(:,3),'-g','linewidth',1)
plot(NEE1(:,4),NEE1(:,5),'-b','linewidth',1)

% plot(NPP_L12(:,4),NPP_L12(:,7),':r','linewidth',3)

set(gca,'ylim',[-8,6]);  %[0 10 20 30 40 50 60 70 80 90]
set(gca,'xlim',[85,960]);  %[0 10 20 30 40 50 60 70 80 90]
clear title xlabel ylabel
set(gca,'xtick',85:90:960,'xticklabel',{'3/25' '6/23' '9/21' '12/20' '3/20' '6/18' '9/16' '12/15' '3/15' '6/13'})
xlabel('Date (mm/dd)')
ylabel('NEE (g C m^{-2} d^{-1})')
legend('OBS','unCPLD','unCPLD-FT','CPLD','Location','North','Orientation','horizontal')
legend boxoff
% text(410,7,'unCPLD: R^2 = 0.74, RMSE= 1.74 g C m^{-2} d^{-1}')
% text(410,6,'CPLD  : R^2 = 0.74, RMSE= 1.77 g C m^{-2} d^{-1}')
% text(0,8,'R^2 = 0.75')
% text(0,7,'RMSE= 0.71 g C m^{-2} d^{-1}')
title('(c)','Position', [85, 6.2, 0])

%% plot the coup/uncoup comparison fig of ecosystem respiration Reco(=GPP+NEE)
load('RECO_UPDATE.mat')
% ind=find(NEE1(:,1)<=0);
% NEE1(ind,1)=NaN;
ind=find(NEE11(:,1)==0);
NEE11(ind,1)=NaN;
ind=find(NEE11(:,4)==0);
NEE11(ind,4)=NaN;
% plot(NEE1(:,4),NEE1(:,1),'om','linewidth',1,'MarkerEdgeColor','m',...
% 'MarkerFaceColor','m')%
figure
plot(NEE11(:,4),NEE11(:,1),'ok','linewidth',1,'MarkerEdgeColor','k',...
'MarkerFaceColor','k','Markersize',4)
hold on
plot(NEE11(:,4),NEE11(:,2),'-r','linewidth',1)
plot(NEE11(:,4),NEE11(:,3),'-g','linewidth',1)
plot(NEE11(:,4),NEE11(:,5),'-b','linewidth',1)

% plot(NPP_L12(:,4),NPP_L12(:,7),':r','linewidth',3)

set(gca,'ylim',[0,15]);  %[0 10 20 30 40 50 60 70 80 90]
set(gca,'xlim',[85,floor(length(QE)/24)+90]);  %[0 10 20 30 40 50 60 70 80 90]
clear title xlabel ylabel
set(gca,'xtick',85:90:(floor(length(QE)/24)+90),'xticklabel',{'3/25' '6/23' '9/21' '12/20' '3/20' '6/18' '9/16' '12/15' '3/15' '6/13'})
xlabel('Date (mm/dd)')
ylabel('R_{eco} (g C m^{-2} d^{-1})')
legend('OBS','unCPLD','unCPLD-FT','CPLD','Location','North','Orientation','horizontal')
legend boxoff
% text(410,7,'unCPLD: R^2 = 0.74, RMSE= 1.74 g C m^{-2} d^{-1}')
% text(410,6,'CPLD  : R^2 = 0.74, RMSE= 1.77 g C m^{-2} d^{-1}')
% text(0,8,'R^2 = 0.75')
% text(0,7,'RMSE= 0.71 g C m^{-2} d^{-1}')
title('(d)','Position', [85, 15.2, 0])

%% Fig 5day average of Rn H LE
%sum up surface flux
%%%% SIMULATED LE
load('MsrH.mat')
load('MsrRn.mat')
load('MsrLE.mat')

NofM=length(QE);
%% unCPLD
    for J=1:floor(NofM/24)
        QES=zeros(1,24);
        QES=QE((J-1)*24+1:(J)*24);
        sumQE(J)=nansum(QES);
    end
        for J=1:floor(NofM/24)
        QES=zeros(1,24);
        QES=Rn((J-1)*24+1:(J)*24);
        sumRn(J)=nansum(QES);
        end
        for J=1:floor(NofM/24)
        QES=zeros(1,24);
        QES=H((J-1)*24+1:(J)*24);
        sumH(J)=nansum(QES);
        end
        
        %% for observations
        H1(H1<-999)=nan;
        for J=1:floor(NofM/24)
        QES=zeros(1,24);
        QES=H1((J-1)*24+1:(J)*24);
        sumH_obs(J)=nansum(QES);
        end
        sumH_obs(sumH_obs(1:251)==0)=nan;
        
        Rn1(Rn1<-999)=nan;
        for J=1:floor(NofM/24)
        QES=zeros(1,24);
        QES=Rn1((J-1)*24+1:(J)*24);
        sumRn_obs(J)=nansum(QES);
        end 
        sumRn_obs(sumRn_obs(1:251)==0)=nan;
        
        QE1(QE1<-999)=nan;
        for J=1:floor(NofM/24)
        QES=zeros(1,24);
        QES=QE1((J-1)*24+1:(J)*24);
        sumQE_obs(J)=nansum(QES);
        end      
        sumQE_obs(sumQE_obs(1:251)==0)=nan;
        %% CPLD
        for J=1:floor(NofM/24)
            QES=zeros(1,24);
            QES=QE((J-1)*24+1:(J)*24);
            sumQE1(J)=nansum(QES);
        end
        for J=1:floor(NofM/24)
            QES=zeros(1,24);
            QES=Rn((J-1)*24+1:(J)*24);
            sumRn1(J)=nansum(QES);
        end
        for J=1:floor(NofM/24)
            QES=zeros(1,24);
            QES=H((J-1)*24+1:(J)*24);
            sumH1(J)=nansum(QES);
        end
         %% unCPLD-FT
        for J=1:floor(NofM/24)
            QES=zeros(1,24);
            QES=QE((J-1)*24+1:(J)*24);
            sumQEFT(J)=nansum(QES);
        end
        for J=1:floor(NofM/24)
            QES=zeros(1,24);
            QES=Rn((J-1)*24+1:(J)*24);
            sumRnFT(J)=nansum(QES);
        end
        for J=1:floor(NofM/24)
            QES=zeros(1,24);
            QES=H((J-1)*24+1:(J)*24);
            sumHFT(J)=nansum(QES);
        end
        
        
        LEH1=QE1+H1;

        for i=1:length(LEH1)
        if abs(Rn1(i))<50 && abs(LEH1(i))>150
            LEH1(i)=nan;
        end
        if abs(Rn1(i))>50 && abs(LEH1(i))/abs(Rn1(i))>2.5
            LEH1(i)=nan;
        end
        end
        
        ind=find(~isnan(LEH1(:,1)));
        % NEE1(ind,1)=NaN;
% ind=find(NEE1(:,4)==0);
% NEE1(ind,4)=NaN;
%% %% energy balance closure plto
%%%%% scatter density plot for LE
figure(6)
ENG_RN=Rn1(ind,1);ENG_LEH=LEH1(ind,1);
hAxes = dscatter(Rn1(ind,1),LEH1(ind,1))
% xlabel(params(1).LongName); ylabel(params(2).LongName);
% hold on
% dscatter(QE,Datam11(:,3),'plottype','contour')
% hold off
% hAxes = dscatter(QE,QE1)
set(gca,'xlim',[-500 1000])
set(gca,'ylim',[-500 1000])
% set(gca,'xtick',-1000:1000:6000,'xticklabel',{-1000 0 1000 2000 3000 4000 5000 6000})
% set(gca,'ytick',-1000:1000:6000,'yticklabel',{-1000 0 1000 2000 3000 4000 5000 6000})
set(gca,'xtick',-500:500:1000,'xticklabel',{-500 0 500 1000})
set(gca,'ytick',-500:500:1000,'yticklabel',{-500 0 500 1000})
p=polyfitZero(Rn1(ind,1),LEH1(ind,1),1)
hold on
h1=refline(1,0) %辅助1:1线
h2=refline(p(1),0) %拟合线获取
set(h1,'LineStyle','--','color','black','linewidth',1.5)
set(h2,'color','blue','linewidth',1.5)
hold off
% set(gca,'xlim',[-150 1000])
% set(gca,'ylim',[-150 1000])
axis on
% mdl = fitlm(Datam11(451:end,1),LEH(451:end))

% p=polyfit(Datam11(1:end,1),Rn(1:end),1)
% hold on
% h1=refline(1,0) %辅助1:1线
% h2=refline(p(1),p(2)) %拟合线获取
% set(h1,'LineStyle','--','color','black','linewidth',1.5)
% set(h2,'color','red','linewidth',1.5)
% hold off
% set(gca,'xlim',[-150 1000])
% set(gca,'ylim',[-150 1000])
% axis on
% mdl = fitlm(Datam11(1:end,1),Rn(1:end))
text(-250,900,'y   = 0.59x')
text(-250,800,'R^2 = 0.85')
text(-470,950,'(a)')

xlabel('Rn (W/m^2)')
ylabel('LE+ H (W/m^2)') % left y-axis 
c=colorbar
colorbar('Ticks',[0,0.2,0.4,0.6,0.8,1])
h_bar(3) = findobj('type','colorbar');
% tick_num = get(h_bar,'ticks');
% tick_str = strsplit(num2str(tick_num*1e8,'%d\n'),'\n');
% set(h_bar,'ticklabels',[0 2 4 6 8]);
h_bar(3).Label.String = 'Frequency';
% set(h_bar,'position',[0.84 0.216 0.03 0.709]);
% set(h_bar,'position',[0.84 0.27 0.03 0.60]);
% change position of '?0^{3}' 
% h_bar.Label.Position=[2.4 0 0];
% h_bar.Label.Rotation = 0 ;
% set(gca,'Position',[.22 .27 .60 .60]);
% title('\rm(c) \it q_L_h')    %title('\rm(c) \it q\rm_L_h')
% set(get(gca,'title'),'Position',[-23 0.8 1.00011])

figure(7)
RNG0=Rn1(ind,1)-G(ind,1);
% hAxes = dscatter(RNG0(:,1),ENG_LEH(:,1))
        for i=1:length(RNG0)
        if abs(RNG0(i))>500 && abs(ENG_LEH(i))<150
            ENG_LEH(i)=nan;
        end
%         if abs(Rn1(i))>50 && abs(ENG_LEH(i))/abs(Rn1(i))>2.5
%             ENG_LEH(i)=nan;
%         end
        end
        
        ind1=find(~isnan(ENG_LEH(:,1)));
ENG_RNG0=RNG0(ind1,1);
ENG_LEH0=ENG_LEH(ind1,1);        
hAxes = dscatter(RNG0(ind1,1),ENG_LEH(ind1,1))
% xlabel(params(1).LongName); ylabel(params(2).LongName);
% hold on
% dscatter(QE,Datam11(:,3),'plottype','contour')
% hold off
% hAxes = dscatter(QE,QE1)
set(gca,'xlim',[-500 1000])
set(gca,'ylim',[-500 1000])
% set(gca,'xtick',-1000:1000:6000,'xticklabel',{-1000 0 1000 2000 3000 4000 5000 6000})
% set(gca,'ytick',-1000:1000:6000,'yticklabel',{-1000 0 1000 2000 3000 4000 5000 6000})
set(gca,'xtick',-500:500:1000,'xticklabel',{-500 0 500 1000})
set(gca,'ytick',-500:500:1000,'yticklabel',{-500 0 500 1000})
p=polyfitZero(RNG0(ind1,1),ENG_LEH(ind1,1),1)
hold on
h1=refline(1,0) %辅助1:1线
h2=refline(p(1),0) %拟合线获取
set(h1,'LineStyle','--','color','black','linewidth',1.5)
set(h2,'color','blue','linewidth',1.5)
hold off
% set(gca,'xlim',[-150 1000])
% set(gca,'ylim',[-150 1000])
axis on
% mdl = fitlm(Datam11(451:end,1),LEH(451:end))

% p=polyfit(Datam11(1:end,1),Rn(1:end),1)
% hold on
% h1=refline(1,0) %辅助1:1线
% h2=refline(p(1),p(2)) %拟合线获取
% set(h1,'LineStyle','--','color','black','linewidth',1.5)
% set(h2,'color','red','linewidth',1.5)
% hold off
% set(gca,'xlim',[-150 1000])
% set(gca,'ylim',[-150 1000])
% axis on
% mdl = fitlm(Datam11(1:end,1),Rn(1:end))
text(-250,900,'y   = 0.70x')
text(-250,800,'R^2 = 0.62')
text(-470,950,'(b)')

xlabel('Rn- G_0 (W/m^2)')
ylabel('LE+ H (W/m^2)') % left y-axis 
c=colorbar
colorbar('Ticks',[0,0.2,0.4,0.6,0.8,1])
h_bar(3) = findobj('type','colorbar');
% tick_num = get(h_bar,'ticks');
% tick_str = strsplit(num2str(tick_num*1e8,'%d\n'),'\n');
% set(h_bar,'ticklabels',[0 2 4 6 8]);
h_bar(3).Label.String = 'Frequency';
% set(h_bar,'position',[0.84 0.216 0.03 0.709]);
% set(h_bar,'position',[0.84 0.27 0.03 0.60]);
% change position of '?0^{3}' 
% h_bar.Label.Position=[2.4 0 0];
% h_bar.Label.Rotation = 0 ;
% set(gca,'Position',[.22 .27 .60 .60]);
% title('\rm(c) \it q_L_h')    %title('\rm(c) \it q\rm_L_h')
% set(get(gca,'title'),'Position',[-23 0.8 1.00011])
% dscatter(QE,Datam11(:,3))
% dscatter
% hAxes = dscatter(QE,QE1)
%%
%% FOR RN
figure(7)
hAxes = dscatter(Rn1(ind,1),LEH1(ind,1))
% xlabel(params(1).LongName); ylabel(params(2).LongName);
% hold on
% dscatter(QE,Datam11(:,3),'plottype','contour')
% hold off
% hAxes = dscatter(QE,QE1)
set(gca,'xlim',[-500 1000])
set(gca,'ylim',[-500 1000])
% set(gca,'xtick',-1000:1000:6000,'xticklabel',{-1000 0 1000 2000 3000 4000 5000 6000})
% set(gca,'ytick',-1000:1000:6000,'yticklabel',{-1000 0 1000 2000 3000 4000 5000 6000})
set(gca,'xtick',-500:500:1000,'xticklabel',{-500 0 500 1000})
set(gca,'ytick',-500:500:1000,'yticklabel',{-500 0 500 1000})
p=polyfitZero(Rn1(ind,1),LEH1(ind,1),1)
hold on
h1=refline(1,0) %辅助1:1线
h2=refline(p(1),0) %拟合线获取
set(h1,'LineStyle','--','color','black','linewidth',1.5)
set(h2,'color','blue','linewidth',1.5)
hold off
% set(gca,'xlim',[-150 1000])
% set(gca,'ylim',[-150 1000])
axis on
mdl = fitlm([sumRn_obs(74:164) sumRn_obs(252:end)]',[sumRn(74:164) sumRn(252:end)]')
text(-500,5500,'R = 0.94')
% text(-50,800,'RMSE= 57.15 W/m^2')
xlabel('Observed Rn (W/m^2)')
ylabel('Simulated Rn (W/m^2)') % left y-axis 
c=colorbar
colorbar('Ticks',[0,0.2,0.4,0.6,0.8,1])
h_bar(3) = findobj('type','colorbar');
h_bar(3).Label.String = 'Frequency';

% LE_AVG_sim=zeros(NofM,24);
% for NoM=MON_LW:MON_HG
%     LE_sim=zeros(1,24);
% 
% %         N=0;d=1;% at 5cm
% %         ST_SUM=0;ST_VAR=0;
% %         for i=(J-1)*24+1:1:24
% %             if Datam(i,2)==NoM
% %                 N=N+1;
% %                 ST_SUM=QE(i)+ST_SUM;
% %                 ST_VAR(N)=QE(i);
% %                 
% %             end
% %         end
% %         LE_sim(J)=ST_SUM/N;
% %         LE_AVG_sim(NoM,J)=LE_sim(J);
% %         LE_VAR(J)=std(ST_VAR);
% %         %         RN_VAR_sim(NoM,J)=RN_VAR(J);
%     end
%     LE_VAR_sim(NoM,:)=LE_VAR(:);
% end
figure(1)
plot(moving(sumH_obs,5),'Color', [200 200 200]./255);
hold on;

plot(moving(sumH,5),'R');
plot(moving(sumH1,5),'G');
plot(moving(sumH11,5),'B');

hold off
mdl = fitlm(sumH_obs(1:end)',sumH(1:end)')
A=[sumH(74:164) sumH(252:end)]';
B=[sumH_obs(74:164) sumH_obs(252:end)]';
BIAS1=sum(A-B)/length(A)

mdl = fitlm(sumH_obs(1:end)',sumH1(1:end)')
A=[sumH1(74:164) sumH1(252:end)]';
B=[sumH_obs(74:164) sumH_obs(252:end)]';
BIAS2=sum(A-B)/length(A)

mdl = fitlm(sumH_obs(1:end)',sumH11(1:end)')
A=[sumH11(74:164) sumH11(252:end)]';
B=[sumH_obs(74:164) sumH_obs(252:end)]';
BIAS3=sum(A-B)/length(A)


xlabel('Date (mm/dd)')
ylabel('H (W/m^2)')
set(gca,'xlim',[1,floor(NofM/24)]);  %[0 10 20 30 40 50 60 70 80 90]
set(gca,'xtick',1:90:floor(NofM/24),'xticklabel',{'3/25' '6/23' '9/21' '12/20' '3/20' '6/18' '9/16' '12/15' '3/15' '6/13'})
set(gca,'ylim',[-500,3500]);  %[0 10 20 30 40 50 60 70 80 90]
set(gca,'ytick',0:1000:3000,'yticklabel',[0 1000 2000 3000])
text(1,3700,'(b)')
legend('OBS','unCPLD','unCPLD-FT','CPLD','Location','North','Orientation','horizontal')
legend boxoff


textstr={'unCPLD: RMSE= 620 '; '                BIAS= 419.61 '};
text(50,2700,textstr,'fontsize',8);
textstr={'unCPLD-FT: RMSE= 602 '; '                      BIAS= 448.03 '};
text(350,2700,textstr,'fontsize',8);
textstr={'CPLD: RMSE= 623 '; '            BIAS= 184.77 '};
text(650,2700,textstr,'fontsize',8);


%% Rn
figure(2)
plot(moving(sumRn_obs,5),'Color', [200 200 200]./255);
hold on;

plot(moving(sumRn,5),'r');
plot(moving(sumRn1,5),'G');
plot(moving(sumRn11,5),'B');
hold off
xlabel('Date (mm/dd)')
ylabel('Rn (W/m^2)')
set(gca,'ylim',[-1000,7000]);  %[0 10 20 30 40 50 60 70 80 90]
set(gca,'xlim',[1,floor(NofM/24)]);  %[0 10 20 30 40 50 60 70 80 90]
set(gca,'xtick',1:90:floor(NofM/24),'xticklabel',{'3/25' '6/23' '9/21' '12/20' '3/20' '6/18' '9/16' '12/15' '3/15' '6/13'})
text(1,7280,'(a)')
legend('OBS','unCPLD','unCPLD-FT','CPLD','Location','North','Orientation','horizontal')
legend boxoff

mdl = fitlm(sumRn_obs(1:end)',sumRn(1:end)')
A=[sumRn(74:164) sumRn(252:end)]';
B=[sumRn_obs(74:164) sumRn_obs(252:end)]';
BIAS1=sum(A-B)/length(A)

mdl = fitlm(sumRn_obs(1:end)',sumRn1(1:end)')
A=[sumRn1(74:164) sumRn1(252:end)]';
B=[sumRn_obs(74:164) sumRn_obs(252:end)]';
BIAS2=sum(A-B)/length(A)

mdl = fitlm(sumRn_obs(1:end)',sumRn11(1:end)')
A=[sumRn11(74:164) sumRn11(252:end)]';
B=[sumRn_obs(74:164) sumRn_obs(252:end)]';
BIAS3=sum(A-B)/length(A)


textstr={'unCPLD: RMSE= 545 '; '                BIAS= 131.30 '};
text(50,5400,textstr,'fontsize',8);
textstr={'unCPLD-FT: RMSE= 531 '; '                      BIAS= 102.71 '};
text(350,5400,textstr,'fontsize',8);
textstr={'CPLD: RMSE= 623 '; '            BIAS= 210.88 '};
text(650,5400,textstr,'fontsize',8);


%% QE
figure(3)
plot(moving(sumQE_obs,5),'Color', [200 200 200]./255);
hold on;

plot(moving(sumQE,5),'r');
plot(moving(sumQE1,5),'G');
plot(moving(sumQE11,5),'B');
hold off
xlabel('Date (mm/dd)')
ylabel('LE (W/m^2)')
set(gca,'xlim',[1,floor(NofM/24)]);  %[0 10 20 30 40 50 60 70 80 90]
set(gca,'xtick',1:90:floor(NofM/24),'xticklabel',{'3/25' '6/23' '9/21' '12/20' '3/20' '6/18' '9/16' '12/15' '3/15' '6/13'})
text(1,4200,'(c)')
set(gca,'ylim',[-500,4000]);  %[0 10 20 30 40 50 60 70 80 90]
set(gca,'ytick',0:1000:4000,'yticklabel',[0 1000 2000 3000 4000])
legend('OBS','unCPLD','unCPLD-FT','CPLD','Location','North','Orientation','horizontal')
legend boxoff


mdl = fitlm(sumQE_obs(1:end)',sumQE(1:end)')
A=[sumQE(74:164) sumQE(252:end)]';
B=[sumQE_obs(74:164) sumQE_obs(252:end)]';
BIAS1=sum(A-B)/length(A)

mdl = fitlm(sumQE_obs(1:end)',sumQE1(1:end)')
A=[sumQE1(74:164) sumQE1(252:end)]';
B=[sumQE_obs(74:164) sumQE_obs(252:end)]';
BIAS2=sum(A-B)/length(A)

mdl = fitlm(sumQE_obs(1:end)',sumQE11(1:end)')
A=[sumQE11(74:164) sumQE11(252:end)]';
B=[sumQE_obs(74:164) sumQE_obs(252:end)]';
BIAS3=sum(A-B)/length(A)

textstr={'unCPLD: RMSE= 397 '; '                BIAS= 49.18 '};
text(50,3100,textstr,'fontsize',8);
textstr={'unCPLD-FT: RMSE= 451 '; '                      BIAS= -7.08 '};
text(350,3100,textstr,'fontsize',8);
textstr={'CPLD: RMSE= 398 '; '            BIAS= 46.31 '};
text(650,3100,textstr,'fontsize',8);




plot(moving(ALB1(:,2),5),'k');
hold on;

plot(moving(ALB1(:,6),5),'r');
plot(moving(ALB1(:,10),5),'b');
hold off
% plot(moving(x,7,'median'),'r');
% plot(moving(x,7,@(x)max(x)),'b');
xlabel('Date (mm/dd)')
ylabel('Rn (W/m^2)')
set(gca,'xlim',[1,270]);  %[0 10 20 30 40 50 60 70 80 90]
set(gca,'xtick',1:60:270,'xticklabel',{'11/20' '1/19' '3/20' '5/19' '7/18'})
text(1,6280,'(a)')
legend('Obs','unCPLD','CPLD','Location','Northwest')
legend boxoff
%%% H
plot(moving(ALB1(:,4),5),'k');
hold on;
plot(moving(ALB1(:,7),5),'r');
plot(moving(ALB1(:,11),5),'b');
hold off
% plot(moving(x,7,'median'),'r');
% plot(moving(x,7,@(x)max(x)),'b');
xlabel('Date (mm/dd)')
ylabel('H (W/m^2)')
set(gca,'xlim',[1,270]);  %[0 10 20 30 40 50 60 70 80 90]
set(gca,'xtick',1:60:270,'xticklabel',{'11/20' '1/19' '3/20' '5/19' '7/18'})
text(1,4200,'(b)')
legend('Obs','unCPLD','CPLD','Location','Northwest')
legend boxoff

%%% LE
plot(moving(ALB1(:,3),5),'k');
hold on;
plot(moving(ALB1(:,8),5),'r');
plot(moving(ALB1(:,12),5),'b');
hold off
% plot(moving(x,7,'median'),'r');
% plot(moving(x,7,@(x)max(x)),'b');
xlabel('Date (mm/dd)')
ylabel('LE (W/m^2)')
set(gca,'ylim',[-500,4000]);  %[0 10 20 30 40 50 60 70 80 90]

set(gca,'xlim',[1,270]);  %[0 10 20 30 40 50 60 70 80 90]
set(gca,'xtick',1:60:270,'xticklabel',{'11/20' '1/19' '3/20' '5/19' '7/18'})
text(1,4200,'(c)')
legend('Obs','unCPLD','CPLD','Location','Northwest')
legend boxoff

%%
%%%%% scatter density plot for LE
% A=sumQE_obs(252:end)';
% n=length(A);
% B=sumQE(252:end)';
% BIAS20=sum(A-B)/length(A)
% %%% RMSE
% RMSE20 = sqrt(mean((A-B).^2))
% [R1,P]=corrcoef(A,B);
figure(4)
hAxes = dscatter([sumQE_obs(74:164) sumQE_obs(252:end)]',[sumQE(74:164) sumQE(252:end)]')
% xlabel(params(1).LongName); ylabel(params(2).LongName);
% hold on
% dscatter(QE,Datam11(:,3),'plottype','contour')
% hold off
% hAxes = dscatter(QE,QE1)
set(gca,'xlim',[-1000 3500])
set(gca,'ylim',[-1000 3500])
set(gca,'xtick',-1000:1000:3000,'xticklabel',{-1000 0 1000 2000 3000})
set(gca,'ytick',-1000:1000:3000,'yticklabel',{-1000 0 1000 2000 3000})

p=polyfit([sumQE_obs(74:164) sumQE_obs(252:end)]',[sumQE(74:164) sumQE(252:end)]',1)
hold on
h1=refline(1,0) %辅助1:1线
h2=refline(p(1),p(2)) %拟合线获取
set(h1,'LineStyle','--','color','black','linewidth',1.5)
set(h2,'color','blue','linewidth',1.5)
hold off
% set(gca,'xlim',[-150 1000])
% set(gca,'ylim',[-150 1000])
axis on
mdl = fitlm(sumQE_obs(1:end)',sumQE(1:end)')
text(-500,3000,'R = 0.89')

xlabel('Observed LE (W/m^2)')
ylabel('Simulated LE (W/m^2)') % left y-axis 
c=colorbar
colorbar('Ticks',[0,0.2,0.4,0.6,0.8,1])
h_bar(3) = findobj('type','colorbar');
% tick_num = get(h_bar,'ticks');
% tick_str = strsplit(num2str(tick_num*1e8,'%d\n'),'\n');
% set(h_bar,'ticklabels',[0 2 4 6 8]);
h_bar(3).Label.String = 'Frequency';
% set(h_bar,'position',[0.84 0.216 0.03 0.709]);
% set(h_bar,'position',[0.84 0.27 0.03 0.60]);
% change position of '?0^{3}' 
% h_bar.Label.Position=[2.4 0 0];
% h_bar.Label.Rotation = 0 ;
% set(gca,'Position',[.22 .27 .60 .60]);
% title('\rm(c) \it q_L_h')    %title('\rm(c) \it q\rm_L_h')
% set(get(gca,'title'),'Position',[-23 0.8 1.00011])
%% unCPLD-FT LE
figure(5)
hAxes = dscatter([sumQE_obs(74:164) sumQE_obs(252:end)]',[sumQE1(74:164) sumQE1(252:end)]')
% xlabel(params(1).LongName); ylabel(params(2).LongName);
% hold on
% dscatter(QE,Datam11(:,3),'plottype','contour')
% hold off
% hAxes = dscatter(QE,QE1)
set(gca,'xlim',[-1000 3500])
set(gca,'ylim',[-1000 3500])
set(gca,'xtick',-1000:1000:3000,'xticklabel',{-1000 0 1000 2000 3000})
set(gca,'ytick',-1000:1000:3000,'yticklabel',{-1000 0 1000 2000 3000})

p=polyfit([sumQE_obs(74:164) sumQE_obs(252:end)]',[sumQE1(74:164) sumQE1(252:end)]',1)
hold on
h1=refline(1,0) %辅助1:1线
h2=refline(p(1),p(2)) %拟合线获取
set(h1,'LineStyle','--','color','black','linewidth',1.5)
set(h2,'color','blue','linewidth',1.5)
hold off
% set(gca,'xlim',[-150 1000])
% set(gca,'ylim',[-150 1000])
axis on
mdl = fitlm([sumQE_obs(74:164) sumQE_obs(252:end)]',[sumQE1(74:164) sumQE1(252:end)]')
text(-500,3000,'R = 0.86')
% text(-50,800,'RMSE= 57.15 W/m^2')

xlabel('Observed LE (W/m^2)')
ylabel('Simulated LE (W/m^2)') % left y-axis 
c=colorbar
colorbar('Ticks',[0,0.2,0.4,0.6,0.8,1])
h_bar(3) = findobj('type','colorbar');
% tick_num = get(h_bar,'ticks');
% tick_str = strsplit(num2str(tick_num*1e8,'%d\n'),'\n');
% set(h_bar,'ticklabels',[0 2 4 6 8]);
h_bar(3).Label.String = 'Frequency';

%% CPLD LE
figure(15)
hAxes = dscatter([sumQE_obs(74:164) sumQE_obs(252:end)]',[sumQE11(74:164) sumQE11(252:end)]')
% xlabel(params(1).LongName); ylabel(params(2).LongName);
% hold on
% dscatter(QE,Datam11(:,3),'plottype','contour')
% hold off
% hAxes = dscatter(QE,QE1)
set(gca,'xlim',[-1000 3500])
set(gca,'ylim',[-1000 3500])
set(gca,'xtick',-1000:1000:3000,'xticklabel',{-1000 0 1000 2000 3000})
set(gca,'ytick',-1000:1000:3000,'yticklabel',{-1000 0 1000 2000 3000})

p=polyfit([sumQE_obs(74:164) sumQE_obs(252:end)]',[sumQE11(74:164) sumQE11(252:end)]',1)
hold on
h1=refline(1,0) %辅助1:1线
h2=refline(p(1),p(2)) %拟合线获取
set(h1,'LineStyle','--','color','black','linewidth',1.5)
set(h2,'color','blue','linewidth',1.5)
hold off
% set(gca,'xlim',[-150 1000])
% set(gca,'ylim',[-150 1000])
axis on
mdl = fitlm([sumQE_obs(74:164) sumQE_obs(252:end)]',[sumQE11(74:164) sumQE11(252:end)]')
text(-500,3000,'R = 0.89')
% text(-50,800,'RMSE= 57.15 W/m^2')

xlabel('Observed LE (W/m^2)')
ylabel('Simulated LE (W/m^2)') % left y-axis 
c=colorbar
colorbar('Ticks',[0,0.2,0.4,0.6,0.8,1])
h_bar(3) = findobj('type','colorbar');
% tick_num = get(h_bar,'ticks');
% tick_str = strsplit(num2str(tick_num*1e8,'%d\n'),'\n');
% set(h_bar,'ticklabels',[0 2 4 6 8]);
h_bar(3).Label.String = 'Frequency';

%% FOR H
figure(6)
hAxes = dscatter([sumH_obs(74:164) sumH_obs(252:end)]',[sumH(74:164) sumH(252:end)]')
% xlabel(params(1).LongName); ylabel(params(2).LongName);
% hold on
% dscatter(QE,Datam11(:,3),'plottype','contour')
% hold off
% hAxes = dscatter(QE,QE1)
set(gca,'xlim',[-1000 3500])
set(gca,'ylim',[-1000 3500])
set(gca,'xtick',-1000:1000:3000,'xticklabel',{-1000 0 1000 2000 3000})
set(gca,'ytick',-1000:1000:3000,'yticklabel',{-1000 0 1000 2000 3000})

p=polyfit([sumH_obs(74:164) sumH_obs(252:end)]',[sumH(74:164) sumH(252:end)]',1)
hold on
h1=refline(1,0) %辅助1:1线
h2=refline(p(1),p(2)) %拟合线获取
set(h1,'LineStyle','--','color','black','linewidth',1.5)
set(h2,'color','blue','linewidth',1.5)
hold off
% set(gca,'xlim',[-150 1000])
% set(gca,'ylim',[-150 1000])
axis on
mdl = fitlm([sumH_obs(74:164) sumH_obs(252:end)]',[sumH(74:164) sumH(252:end)]')
text(-500,3000,'R = 0.49')
% text(-50,800,'RMSE= 57.15 W/m^2')
xlabel('Observed H (W/m^2)')
ylabel('Simulated H (W/m^2)') % left y-axis 
c=colorbar
colorbar('Ticks',[0,0.2,0.4,0.6,0.8,1])
h_bar(3) = findobj('type','colorbar');
h_bar(3).Label.String = 'Frequency';

figure(7)
hAxes = dscatter([sumH_obs(74:164) sumH_obs(252:end)]',[sumH1(74:164) sumH1(252:end)]')
% xlabel(params(1).LongName); ylabel(params(2).LongName);
% hold on
% dscatter(QE,Datam11(:,3),'plottype','contour')
% hold off
% hAxes = dscatter(QE,QE1)
set(gca,'xlim',[-1000 3500])
set(gca,'ylim',[-1000 3500])
set(gca,'xtick',-1000:1000:3000,'xticklabel',{-1000 0 1000 2000 3000})
set(gca,'ytick',-1000:1000:3000,'yticklabel',{-1000 0 1000 2000 3000})

p=polyfit([sumH_obs(74:164) sumH_obs(252:end)]',[sumH1(74:164) sumH1(252:end)]',1)
hold on
h1=refline(1,0) %辅助1:1线
h2=refline(p(1),p(2)) %拟合线获取
set(h1,'LineStyle','--','color','black','linewidth',1.5)
set(h2,'color','blue','linewidth',1.5)
hold off
% set(gca,'xlim',[-150 1000])
% set(gca,'ylim',[-150 1000])
axis on
mdl = fitlm([sumH_obs(74:164) sumH_obs(252:end)]',[sumH1(74:164) sumH1(252:end)]')
text(-500,3000,'R = 0.50')
% text(-50,800,'RMSE= 57.15 W/m^2')
xlabel('Observed H (W/m^2)')
ylabel('Simulated H (W/m^2)') % left y-axis 
c=colorbar
colorbar('Ticks',[0,0.2,0.4,0.6,0.8,1])
h_bar(3) = findobj('type','colorbar');
h_bar(3).Label.String = 'Frequency';

figure(17)
hAxes = dscatter([sumH_obs(74:164) sumH_obs(252:end)]',[sumH11(74:164) sumH11(252:end)]')
% xlabel(params(1).LongName); ylabel(params(2).LongName);
% hold on
% dscatter(QE,Datam11(:,3),'plottype','contour')
% hold off
% hAxes = dscatter(QE,QE1)
set(gca,'xlim',[-1000 3500])
set(gca,'ylim',[-1000 3500])
set(gca,'xtick',-1000:1000:3000,'xticklabel',{-1000 0 1000 2000 3000})
set(gca,'ytick',-1000:1000:3000,'yticklabel',{-1000 0 1000 2000 3000})

p=polyfit([sumH_obs(74:164) sumH_obs(252:end)]',[sumH11(74:164) sumH11(252:end)]',1)
hold on
h1=refline(1,0) %辅助1:1线
h2=refline(p(1),p(2)) %拟合线获取
set(h1,'LineStyle','--','color','black','linewidth',1.5)
set(h2,'color','blue','linewidth',1.5)
hold off
% set(gca,'xlim',[-150 1000])
% set(gca,'ylim',[-150 1000])
axis on
mdl = fitlm([sumH_obs(74:164) sumH_obs(252:end)]',[sumH11(74:164) sumH11(252:end)]')
text(-500,3000,'R = 0.49')
% text(-50,800,'RMSE= 57.15 W/m^2')
xlabel('Observed H (W/m^2)')
ylabel('Simulated H (W/m^2)') % left y-axis 
c=colorbar
colorbar('Ticks',[0,0.2,0.4,0.6,0.8,1])
h_bar(3) = findobj('type','colorbar');
h_bar(3).Label.String = 'Frequency';
%% FOR RN
figure(7)
hAxes = dscatter([sumRn_obs(74:164) sumRn_obs(252:end)]',[sumRn(74:164) sumRn(252:end)]')
% xlabel(params(1).LongName); ylabel(params(2).LongName);
% hold on
% dscatter(QE,Datam11(:,3),'plottype','contour')
% hold off
% hAxes = dscatter(QE,QE1)
set(gca,'xlim',[-1000 6000])
set(gca,'ylim',[-1000 6000])
% set(gca,'xtick',-1000:1000:6000,'xticklabel',{-1000 0 1000 2000 3000 4000 5000 6000})
% set(gca,'ytick',-1000:1000:6000,'yticklabel',{-1000 0 1000 2000 3000 4000 5000 6000})
set(gca,'xtick',0:2000:6000,'xticklabel',{0 2000 4000 6000})
set(gca,'ytick',0:2000:6000,'yticklabel',{0 2000 4000 6000})
p=polyfit([sumRn_obs(74:164) sumRn_obs(252:end)]',[sumRn(74:164) sumRn(252:end)]',1)
hold on
h1=refline(1,0) %辅助1:1线
h2=refline(p(1),p(2)) %拟合线获取
set(h1,'LineStyle','--','color','black','linewidth',1.5)
set(h2,'color','blue','linewidth',1.5)
hold off
% set(gca,'xlim',[-150 1000])
% set(gca,'ylim',[-150 1000])
axis on
mdl = fitlm([sumRn_obs(74:164) sumRn_obs(252:end)]',[sumRn(74:164) sumRn(252:end)]')
text(-500,5500,'R = 0.94')
% text(-50,800,'RMSE= 57.15 W/m^2')
xlabel('Observed Rn (W/m^2)')
ylabel('Simulated Rn (W/m^2)') % left y-axis 
c=colorbar
colorbar('Ticks',[0,0.2,0.4,0.6,0.8,1])
h_bar(3) = findobj('type','colorbar');
h_bar(3).Label.String = 'Frequency';

figure(8)
hAxes = dscatter([sumRn_obs(74:164) sumRn_obs(252:end)]',[sumRn1(74:164) sumRn1(252:end)]')
% xlabel(params(1).LongName); ylabel(params(2).LongName);
% hold on
% dscatter(QE,Datam11(:,3),'plottype','contour')
% hold off
% hAxes = dscatter(QE,QE1)
set(gca,'xlim',[-1000 6000])
set(gca,'ylim',[-1000 6000])
% set(gca,'xtick',-1000:1000:6000,'xticklabel',{-1000 0 1000 2000 3000 4000 5000 6000})
% set(gca,'ytick',-1000:1000:6000,'yticklabel',{-1000 0 1000 2000 3000 4000 5000 6000})
set(gca,'xtick',0:2000:6000,'xticklabel',{0 2000 4000 6000})
set(gca,'ytick',0:2000:6000,'yticklabel',{0 2000 4000 6000})
p=polyfit([sumRn_obs(74:164) sumRn_obs(252:end)]',[sumRn1(74:164) sumRn1(252:end)]',1)
hold on
h1=refline(1,0) %辅助1:1线
h2=refline(p(1),p(2)) %拟合线获取
set(h1,'LineStyle','--','color','black','linewidth',1.5)
set(h2,'color','blue','linewidth',1.5)
hold off
% set(gca,'xlim',[-150 1000])
% set(gca,'ylim',[-150 1000])
axis on
mdl = fitlm([sumRn_obs(74:164) sumRn_obs(252:end)]',[sumRn1(74:164) sumRn1(252:end)]')
text(-500,5500,'R = 0.94')
% text(-50,800,'RMSE= 57.15 W/m^2')
xlabel('Observed Rn (W/m^2)')
ylabel('Simulated Rn (W/m^2)') % left y-axis 
c=colorbar
colorbar('Ticks',[0,0.2,0.4,0.6,0.8,1])
h_bar(3) = findobj('type','colorbar');
h_bar(3).Label.String = 'Frequency';


figure(8)
hAxes = dscatter([sumRn_obs(74:164) sumRn_obs(252:end)]',[sumRn11(74:164) sumRn11(252:end)]')
% xlabel(params(1).LongName); ylabel(params(2).LongName);
% hold on
% dscatter(QE,Datam11(:,3),'plottype','contour')
% hold off
% hAxes = dscatter(QE,QE1)
set(gca,'xlim',[-1000 6000])
set(gca,'ylim',[-1000 6000])
% set(gca,'xtick',-1000:1000:6000,'xticklabel',{-1000 0 1000 2000 3000 4000 5000 6000})
% set(gca,'ytick',-1000:1000:6000,'yticklabel',{-1000 0 1000 2000 3000 4000 5000 6000})
set(gca,'xtick',0:2000:6000,'xticklabel',{0 2000 4000 6000})
set(gca,'ytick',0:2000:6000,'yticklabel',{0 2000 4000 6000})
p=polyfit([sumRn_obs(74:164) sumRn_obs(252:end)]',[sumRn11(74:164) sumRn11(252:end)]',1)
hold on
h1=refline(1,0) %辅助1:1线
h2=refline(p(1),p(2)) %拟合线获取
set(h1,'LineStyle','--','color','black','linewidth',1.5)
set(h2,'color','blue','linewidth',1.5)
hold off
% set(gca,'xlim',[-150 1000])
% set(gca,'ylim',[-150 1000])
axis on
mdl = fitlm([sumRn_obs(74:164) sumRn_obs(252:end)]',[sumRn11(74:164) sumRn11(252:end)]')
text(-500,5500,'R = 0.92')
% text(-50,800,'RMSE= 57.15 W/m^2')
xlabel('Observed Rn (W/m^2)')
ylabel('Simulated Rn (W/m^2)') % left y-axis 
c=colorbar
colorbar('Ticks',[0,0.2,0.4,0.6,0.8,1])
h_bar(3) = findobj('type','colorbar');
h_bar(3).Label.String = 'Frequency';

%% ice-flux plotting
x=[1:1:length(QE)];
% y2=[10 20 50 100 150 200 300 400 500 600 800 1000 1250 1500 1750 2000 2500 3000]./10;
% y3(1:1:length(y2))=y2(length(y2):-1:1);
y3=[18:-1:1];
TTICE=TSim_Theta_I';
QTTICE(:,1:1:18)=TTICE(:,18:-1:1);
QLLPT=QLHHT'.*36000;
QLLPT(abs(QLLPT)>2.5)=0;
QLLPT(:,17:18)=nan;

Qsave(:,1:1:17)=qsave(:,17:-1:1);
Qsave(:,17:18)=nan;
Qsave(Qsave>2.5)=0;
Qsave=-Qsave;

subplot(3,1,1)
%%%SOIL ICE CONTENT

h=imagesc(x,y3,QTTICE');
% set(gca,'Ydir','normal');
% h=imagesc(QTTICE')
set(h,'alphadata',~isnan(QTTICE'))
set(gca,'xtick',1:90*24:length(QE))
set(gca,'XTickLabel',[]);
% xlabel('Date (mm/dd)'),'xticklabel',{'11/20' '1/19' '3/20' '5/19' '7/18'}
set(gca,'ytick',[1 5 10 15],'yticklabel',[1 10 50 150])

ylabel('Soil depth (cm)')
title('(a) Soil ice content')%,'Position', [0, 0.5, 0]
hold on
set(gca,'ylim',[1,18]);  %[0 10 20 30 40 50 60 70 80 90]
set(gca,'xlim',[0,length(QE)]);  %[0 10 20 30 40 50 60 70 80 90]

ylim=get(gca,'Ylim'); 
for i=1:floor(length(QE)/24/90)
plot([i*90*24,i*90*24],ylim,'k--'); % ??x=8,9,10,11???
end
c=colorbar
colorbar('Ticks',[0,0.05,0.1,0.15])
%%%%%
subplot(3,1,3)
%%%%% COUPLED
% figure
% [c1,h1]=contourf(x,y2,QLLPT');%,'Color','k','ShowText','on'
% set(h1,'ShowText','on','TextList',[0,0])
% set(gca,'Ydir','reverse');
h=imagesc(x,y3,QLLPT');
set(h,'alphadata',~isnan(QLLPT'))

colormap(jet);
set(gca,'xtick',1:90*24:length(QE),'xticklabel',{'3/25' '6/23' '9/21' '12/20' '3/20' '6/18' '9/16' '12/15' '3/15' '6/13'})
% xlabel('Date (mm/dd)'),'xticklabel',{'11/20' '1/19' '3/20' '5/19' '7/18'}
% set(gca,'XTickLabel',[]);
% xlabel('Date (mm/dd)'),'xticklabel',{'11/20' '1/19' '3/20' '5/19' '7/18'}
set(gca,'ytick',[1 5 10 15],'yticklabel',[1 10 50 150])
xlabel('Date (mm/dd)')
ylabel('Soil depth (cm)')
title('(c) CPLD')%,'Position', [0, 0.5, 0]
hold on
set(gca,'ylim',[1,18]);  %[0 10 20 30 40 50 60 70 80 90]
set(gca,'xlim',[0,length(QE)]);  %[0 10 20 30 40 50 60 70 80 90]

ylim=get(gca,'Ylim'); 
for i=1:floor(length(QE)/24/90)
plot([i*90*24,i*90*24],ylim,'k--'); % ??x=8,9,10,11???
end
c=colorbar
colorbar('Ticks',[-2,-1,0,1,2])
lim = caxis
caxis([-1.8621 2.05])
%%% uncoupled T&C
subplot(3,1,2)

h=imagesc(x,y3,Qsave');
set(h,'alphadata',~isnan(Qsave'))

% set(gca,'Ydir','normal');
set(gca,'xtick',1:90*24:length(QE),'xticklabel',{'3/25' '6/23' '9/21' '12/20' '3/20' '6/18' '9/16' '12/15' '3/15' '6/13'})
set(gca,'XTickLabel',[]);
% xlabel('Date (mm/dd)')
title('(b) unCPLD')%,'Position', [0, 0.5, 0]
% set(gca,'XTickLabel',[]);
% xlabel('Date (mm/dd)'),'xticklabel',{'11/20' '1/19' '3/20' '5/19' '7/18'}
set(gca,'ytick',[1 5 10 15],'yticklabel',[1 10 50 150])

ylabel('Soil depth (cm)')
title('(b) unCPLD')%,'Position', [0, 0.5, 0]
hold on
set(gca,'ylim',[1,18]);  %[0 10 20 30 40 50 60 70 80 90]
set(gca,'xlim',[0,length(QE)]);  %[0 10 20 30 40 50 60 70 80 90]
lim = caxis
caxis([-1.8621 2.05])

ylim=get(gca,'Ylim'); 
for i=1:floor(length(QE)/24/90)
plot([i*90*24,i*90*24],ylim,'k--'); % ??x=8,9,10,11???
end
h=colorbar
colorbar('Ticks',[-2,-1,0,1,2])
% % h = colorbar;
% ylabel(h, 'mm/h')
h_bar = findobj('type','colorbar');
h_bar(1).Label.String ='Flux (mm/h)';
h_bar(2).Label.String ='Flux (mm/h)';
h_bar(3).Label.String ='Ice content (-)';

% % c=colorbar
% % colorbar('Ticks',[0,2e-8,4e-8,6e-8,8e-8])
% 
% tick_num = get(h_bar,'ticks');
% tick_str = strsplit(num2str(tick_num*1e8,'%d\n'),'\n');
% set(h_bar,'ticklabels',[0 2 4 6 8]);
% % set(h_bar,'position',[0.84 0.216 0.03 0.709]);
% set(h_bar,'position',[0.84 0.27 0.03 0.60]);
% % change position of '?0^{3}' 
% h_bar.Label.Position=[2.4 0 0];
% h_bar.Label.Rotation = 0 ;


%% plot the coup/uncoup comparison fig of GPP
load('NEE_COMPUPDATE.mat')
% ind=find(NEE1(:,1)<=0);
% NEE1(ind,1)=NaN;
ind=find(NEE1(:,1)==0);
NEE1(ind,1)=NaN;
ind=find(NEE1(:,4)==0);
NEE1(ind,4)=NaN;
% plot(NEE1(:,4),NEE1(:,1),'om','linewidth',1,'MarkerEdgeColor','m',...
% 'MarkerFaceColor','m')%

plot(NEE1(:,4),NEE1(:,1),'*-K','linewidth',1,'Markersize',8)
hold on
plot(NEE1(:,4),NEE1(:,2),'-b','linewidth',1)
plot(NEE1(:,4),NEE1(:,3),'-r','linewidth',1)
% plot(NPP_L12(:,4),NPP_L12(:,7),':r','linewidth',3)

set(gca,'ylim',[-10,5]);  %[0 10 20 30 40 50 60 70 80 90]
set(gca,'xlim',[85,960]);  %[0 10 20 30 40 50 60 70 80 90]
clear title xlabel ylabel
set(gca,'xtick',85:90:960,'xticklabel',{'3/25' '6/23' '9/21' '12/20' '3/20' '6/18' '9/16' '12/15' '3/15' '6/13'})
xlabel('Date (mm/dd)')
ylabel('NEE (g C m^{-2} d^{-1})')
legend('OBS','unCPLD','CPLD','Location','Southwest')
legend boxoff
% text(410,7,'unCPLD: R^2 = 0.74, RMSE= 1.74 g C m^{-2} d^{-1}')
% text(410,6,'CPLD  : R^2 = 0.74, RMSE= 1.77 g C m^{-2} d^{-1}')
% text(0,8,'R^2 = 0.75')
% text(0,7,'RMSE= 0.71 g C m^{-2} d^{-1}')
% title('(a)','Position', [324, 18.2, 0])



set(gca,'fontsize',22);
% set(gca,'yticklabel',[0.1 0.5 1 2]);  %,'xticklabel',[0 0.5 0.8 0.9 1.5 1.8 2][0 10 20 30 40 50 60 70 80 90]

% % contourf(XX,YY,VPtop_QL_HH')
hold on
% [c,h1]=contour(XX,YY,VP_SimQL_HH,[-2e-8 -1e-8 0 1e-8 2e-8],'LineWidth',2,'Color','k');
% set(h1,'ShowText','on','LevelList',[-2e-8 -1e-8 0 1e-8 2e-8]);
% figure
% d =5;dd =1;%%%??d????????dd?????????????,'LineWidth',1
% h=quiver(XX(1:end,1:d:end),YY(1:end,1:d:end),FX(1:end,1:d:end)./dd,FY(1:end,1:d:end)./dd,0.5)
% 
% % h=quiver(XX,YY,FX,FY,0.8)
% set(h,'maxheadsize',5,'color','r');
% ax=gca;
% set(gca,'Ydir','reverse');

% set(h,'maxheadsize',5.2);
hold off
ylabel('Soil depth (cm)')
% xlabel('Time (30 min)')
% set(gca,'fontsize',14);
% set(gca,'xtick',0:10:5067/48);
% set(gca,'xtick',410:529,'xticklabel',0:0.5:72/48); 
% set(gca,'xtick',410:24:529,'xticklabel',0:0.5:119/24);
% xlabel('Day of year')

set(gca,'xtick',0:48:240,'xticklabel',7:1:240/48+7); 
% set(gca,'xcolor','none'); 
xlabel('Day after Dec. 1, 2015','Color','r')

set(gca,'ytick',[1 5 6],'yticklabel',[0.1 1 2]);  %,'xticklabel',[0 0.5 0.8 0.9 1.5 1.8 2][0 10 20 30 40 50 60 70 80 90]
set(gca,'ylim',[1,6]);  %[0 10 20 30 40 50 60 70 80 90]
ylabel('Soil depth (cm)')
hold on
ylim=get(gca,'Ylim'); 
for i=1:4
plot([i*48,i*48],ylim,'k--'); % ??x=8,9,10,11???
end
% set(gca,'position',[0.15 0.216 0.67 0.709]);

% c=colorbar
% colorbar('Ticks',[0,2e-8,4e-8,6e-8,8e-8])
h_bar = findobj('type','colorbar');
% tick_num = get(h_bar,'ticks');
% tick_str = strsplit(num2str(tick_num*1e8,'%d\n'),'\n');
% set(h_bar,'ticklabels',[0 2 4 6 8]);
% h_bar.Label.String = '?0^{-8}';
% set(h_bar,'position',[0.84 0.216 0.03 0.709]);
set(h_bar,'position',[0.84 0.27 0.03 0.60]);
% change position of '?0^{3}' 
% h_bar.Label.Position=[2.4 0 0];
% h_bar.Label.Rotation = 0 ;
set(gca,'Position',[.22 .27 .60 .60]);
title('\rm(c) \it q_L_h')    %title('\rm(c) \it q\rm_L_h')
set(get(gca,'title'),'Position',[-23 0.8 1.00011])



%% for the UNcoupled version 
%%%%% scatter density plot for LE
figure(4)
hAxes = dscatter(Datam11(:,3),QE)
% xlabel(params(1).LongName); ylabel(params(2).LongName);
% hold on
% dscatter(QE,Datam11(:,3),'plottype','contour')
% hold off
p=polyfit(Datam11(451:end,3),QE(451:end),1)
hold on
h1=refline(1,0) %辅助1:1线
h2=refline(p(1),p(2)) %拟合线获取
set(h1,'LineStyle','--','color','black','linewidth',1.5)
set(h2,'color','red','linewidth',1.5)
hold off
% hAxes = dscatter(QE,QE1)
set(gca,'xlim',[-100 900])
set(gca,'ylim',[-100 900])
mdl = fitlm(Datam11(451:end,3),QE(451:end))
set(gca,'xtick',0:200:900,'xticklabel',{0 200 400 600 800})
set(gca,'ytick',0:200:900,'yticklabel',{0 200 400 600 800})
text(0,800,'R^2 = 0.73')
text(0,700,'RMSE= 51 W/m^2')

xlabel('Observed LE (W/m^2)')
ylabel('Simulated LE (W/m^2)') % left y-axis 
c=colorbar;
colorbar('Ticks',[0,0.2,0.4,0.6,0.8,1])
h_bar(1) = findobj('type','colorbar');
% tick_num = get(h_bar,'ticks');
% tick_str = strsplit(num2str(tick_num*1e8,'%d\n'),'\n');
% set(h_bar,'ticklabels',[0 2 4 6 8]);
h_bar(1).Label.String = 'Frequency';
% set(h_bar,'position',[0.84 0.216 0.03 0.709]);
% set(h_bar,'position',[0.84 0.27 0.03 0.60]);
% change position of '?0^{3}' 
% h_bar.Label.Position=[2.4 0 0];
% h_bar.Label.Rotation = 0 ;
% set(gca,'Position',[.22 .27 .60 .60]);
% title('\rm(c) \it q_L_h')    %title('\rm(c) \it q\rm_L_h')
% set(get(gca,'title'),'Position',[-23 0.8 1.00011])
% clear h_bar

%%%%% scatter density plot for LE
figure(5)
hAxes = dscatter(Datam11(451:end,2),H(451:end))
% xlabel(params(1).LongName); ylabel(params(2).LongName);
% hold on
% dscatter(QE,Datam11(:,3),'plottype','contour')
% hold off
% hAxes = dscatter(QE,QE1)
set(gca,'xlim',[-150 800])
set(gca,'ylim',[-150 800])
set(gca,'xtick',0:200:800,'xticklabel',{0 200 400 600 800})
set(gca,'ytick',0:200:800,'yticklabel',{0 200 400 600 800})

p=polyfit(Datam11(451:end,2),H(451:end),1)
hold on
h1=refline(1,0) %辅助1:1线
h2=refline(p(1),p(2)) %拟合线获取
set(h1,'LineStyle','--','color','black','linewidth',1.5)
set(h2,'color','red','linewidth',1.5)
hold off
set(gca,'xlim',[-150 800])
set(gca,'ylim',[-150 800])
axis on
mdl = fitlm(Datam11(451:end,2),H(451:end))
text(-50,700,'R^2 = 0.66')
text(-50,600,'RMSE= 57.80 W/m^2')

xlabel('Observed H (W/m^2)')
ylabel('Simulated H (W/m^2)') % left y-axis 
c=colorbar
colorbar('Ticks',[0,0.2,0.4,0.6,0.8,1])
h_bar(2) = findobj('type','colorbar');
% tick_num = get(h_bar,'ticks');
% tick_str = strsplit(num2str(tick_num*1e8,'%d\n'),'\n');
% set(h_bar,'ticklabels',[0 2 4 6 8]);
h_bar(2).Label.String = 'Frequency';
% set(h_bar,'position',[0.84 0.216 0.03 0.709]);
% set(h_bar,'position',[0.84 0.27 0.03 0.60]);
% change position of '?0^{3}' 
% h_bar.Label.Position=[2.4 0 0];
% h_bar.Label.Rotation = 0 ;
% set(gca,'Position',[.22 .27 .60 .60]);
% title('\rm(c) \it q_L_h')    %title('\rm(c) \it q\rm_L_h')
% set(get(gca,'title'),'Position',[-23 0.8 1.00011])



