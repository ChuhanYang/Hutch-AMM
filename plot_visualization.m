% This code is for generating and tuning visualization only

% Visualization: Estimation error comparison
% load('workspace_saved/ExpDec_h10_ErrEst')
% load('workspace_saved/LinInc_h10_ErrEst')
% load('workspace_saved/Drift_h10_ErrEst')
% load('workspace_saved/Uniform_h10_ErrEst')

% f = figure();
% result_tbl = [mean(result_tbl_h1);mean(result_tbl_h3);mean(result_tbl_h5); mean(result_tbl_h10);mean(result_tbl_h20);mean(result_tbl_h50)];
% plot(s_list,result_tbl,'o-','LineWidth', 1,'MarkerSize',10)
% %ylim([0,0.18])
% xlabel('AMM sample number','FontSize',30);
% ylabel('mean relative error','FontSize',30);
% title('Block Size = 20, Estimation error comparison','fontweight','bold','FontSize',36,'interpreter','latex')
% hold on;
% plot(s_list,mean(result_tbl_opt),'*-','LineWidth', 1,'MarkerSize',10)
% plot(s_list,mean(result_tbl_uni),'*-','LineWidth', 1,'MarkerSize',10)
% legend('Hutch AMM h = 1','Hutch AMM h = 3','Hutch AMM h = 5','Hutch AMM h = 10','Hutch AMM h = 20','Hutch AMM h = 50','Optimal AMM','Uniform Sampling','FontSize',20,'interpreter','latex')
% hold off;


f = figure();
result_tbl = [mean(result_tbl_h1); mean(result_tbl_h10);mean(result_tbl_h50)];
plot(s_list,result_tbl,'o-','LineWidth', 1,'MarkerSize',10)
%ylim([0,0.2])
xlabel('AMM sample number','FontSize',30);
ylabel('mean relative error','FontSize',30);
title('Uniform, Estimation error comparison','fontweight','bold','FontSize',36,'interpreter','latex')
hold on;
plot(s_list,mean(result_tbl_opt),'*-','LineWidth', 1,'MarkerSize',10,'Color','#4B0092')
plot(s_list,mean(result_tbl_uni),'*-','LineWidth', 1,'MarkerSize',10,'Color','#40B0A6')
set(gca, 'YScale', 'log')
legend('Hutch AMM h = 1','Hutch AMM h = 10','Hutch AMM h = 50','Optimal AMM','Uniform Sampling','FontSize',20,'interpreter','latex')
hold off;

% Requires R2020a or later
%exportgraphics(f,'Error1.pdf','Resolution',600) 

% Visualization: Estimation error comparison VS time measurements
% load('workspace_saved/error_vs_time')
% load('workspace_saved/error_vs_time_bsz100')

f = figure();
hold on
for i = 1:1:repeat_time
    scatter(time_tbl1,result_tbl_opt(i,:),80,'filled','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
    scatter(time_tbl3,result_tbl_uni(i,:),80,'filled','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
    scatter(time_tbl2,result_tbl_h10(i,:),80,'filled','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
end
[t1,I1] = sort(time_tbl1,'ascend');
[t2,I2] = sort(time_tbl2,'ascend');
[t3,I3] = sort(time_tbl3,'ascend');
plot(t1,mean(result_tbl_opt(:,I1)),'Color','b','LineWidth', 1)
plot(t3,mean(result_tbl_uni(:,I3)),'Color','r','LineWidth', 1)
plot(t2,mean(result_tbl_h10(:,I2)),'Color','g','LineWidth', 1)
set(gca, 'YScale', 'log')
xlabel('time (s)','FontSize',30);
ylabel('relative error','FontSize',30);
title('Increasing sampling number c: reconstruction performance vs. computation time','fontweight','bold','FontSize',36,'interpreter','latex')
legend('Block AMM','Uniform Sampling','Hutch h=10','FontSize',20,'interpreter','latex','Location','northeast')
hold off

% Visualization: Estimation error comparison VS time measurements (3rd version)
load('workspace_saved/error_vs_time_s500_bsz100_500x1m')
t_tbl_AMM = t_tbl_AMM(1:5,:);
result_tbl_AMM = result_tbl_AMM(1:5,:);

f = figure();
hold on
for i = 1:1:repeat_time
    scatter(t_tbl_opt(i,:),result_tbl_opt(i,:),80,'filled','MarkerFaceColor',[255 194 10]/255,'MarkerEdgeColor',[255 194 10]/255,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.5)
    scatter(t_tbl_uni(i,:),result_tbl_uni(i,:),80,'filled','MarkerFaceColor',[26 255 26]/255,'MarkerEdgeColor',[26 255 26]/255,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.5)
    scatter(t_tbl_h10(i,:),result_tbl_h10(i,:),80,'filled','MarkerFaceColor',[12 123 220]/255,'MarkerEdgeColor',[12 123 220]/255,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.5)
    scatter(t_tbl_AMM(i,:),result_tbl_AMM(i,:),80,'filled','MarkerFaceColor',[75 0 146]/255,'MarkerEdgeColor',[75 0 146]/255,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.5)
end
[t1,I1] = sort(mean(t_tbl_opt),'ascend');
[t2,I2] = sort(mean(t_tbl_h10),'ascend');
[t3,I3] = sort(mean(t_tbl_uni),'ascend');
[t4,I4] = sort(mean(t_tbl_AMM),'ascend');
plot(t1,mean(result_tbl_opt(:,I1)),'Color',[255 194 10]/255,'LineWidth', 3)
plot(t3,mean(result_tbl_uni(:,I3)),'Color',[26 255 26]/255,'LineWidth', 3)
plot(t2,mean(result_tbl_h10(:,I2)),'Color',[12 123 220]/255,'LineWidth', 3)
plot(t4,mean(result_tbl_AMM(:,I4)),'Color',[75 0 146]/255,'LineWidth', 3)
xline(std_time,'--','Standard Matrix Multiplication','LineWidth',1.5,'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center'); %11.3601202334
set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')
xlim([0,4]);
ylim([0,1.1]);
xlabel('time (s)','FontSize',30);
ylabel('relative error','FontSize',30);
title('Increasing sampling number c: reconstruction performance vs. computation time','fontweight','bold','FontSize',36,'interpreter','latex')
legend('Block AMM','Uniform Sampling','Hutch h=10','Standard AMM','FontSize',20,'interpreter','latex','Location','northeast')
hold off

%Visualization: Bubble chart
t = tiledlayout(1,1);
nexttile
bubblechart(mean(t_tbl_opt),mean(result_tbl_opt),std(t_tbl_opt),'#FFC20A','MarkerFaceAlpha',0.2)
hold on
bubblechart(mean(t_tbl_uni),mean(result_tbl_uni),std(t_tbl_uni),'#1AFF1A','MarkerFaceAlpha',0.2)
bubblechart(mean(t_tbl_h10),mean(result_tbl_h10),std(t_tbl_h10),'#0C7BDC','MarkerFaceAlpha',0.2)
bubblechart(mean(t_tbl_AMM),mean(result_tbl_AMM),std(t_tbl_AMM),'#4B0092','MarkerFaceAlpha',0.2)
bubblelegend('time standard deviation','FontSize',14,'Location','southeast');
% [t1,I1] = sort(mean(t_tbl_opt),'ascend');
% [t2,I2] = sort(mean(t_tbl_h10),'ascend');
% [t3,I3] = sort(mean(t_tbl_uni),'ascend');
% [t4,I4] = sort(mean(t_tbl_AMM),'ascend');
% plot(t1,mean(result_tbl_opt(:,I1)),'Color',[255 194 10]/255,'LineWidth', 3)
% plot(t3,mean(result_tbl_uni(:,I3)),'Color',[26 255 26]/255,'LineWidth', 3)
% plot(t2,mean(result_tbl_h10(:,I2)),'Color',[12 123 220]/255,'LineWidth', 3)
% plot(t4,mean(result_tbl_AMM(:,I4)),'Color',[75 0 146]/255,'LineWidth', 3)
xline(std_time,'--','Standard Matrix Multiplication','LineWidth',1.5,'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center');
set(gca, 'YScale', 'log')
xlabel('time (s)','FontSize',30);
ylabel('relative error','FontSize',30);
title('Bubble chart','fontweight','bold','FontSize',36,'interpreter','latex')
legend('Block AMM','Uniform Sampling','Hutch h=10','Standard AMM','FontSize',20,'interpreter','latex','Location','northeast')
hold off

%Visualization: Bubble chart 2
t = tiledlayout(1,1);
nexttile
hold on
for i = 1:1:repeat_time
    scatter(t_tbl_opt(i,:),result_tbl_opt(i,:),80,'filled','+','MarkerFaceColor',[255 194 10]/255,'MarkerEdgeColor',[255 194 10]/255,'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',1)
    scatter(t_tbl_uni(i,:),result_tbl_uni(i,:),80,'filled','pentagram','MarkerFaceColor',[26 255 26]/255,'MarkerEdgeColor',[26 255 26]/255,'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',1)
    scatter(t_tbl_h10(i,:),result_tbl_h10(i,:),80,'filled','square','MarkerFaceColor',[12 123 220]/255,'MarkerEdgeColor',[12 123 220]/255,'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',1)
    scatter(t_tbl_AMM(i,:),result_tbl_AMM(i,:),80,'filled','*','MarkerFaceColor',[75 0 146]/255,'MarkerEdgeColor',[75 0 146]/255,'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',1)
end
bubblechart(mean(t_tbl_opt),mean(result_tbl_opt),var(t_tbl_opt),'#FFC20A','MarkerFaceAlpha',0.2)
bubblechart(mean(t_tbl_uni),mean(result_tbl_uni),var(t_tbl_uni),'#1AFF1A','MarkerFaceAlpha',0.2)
bubblechart(mean(t_tbl_h10),mean(result_tbl_h10),var(t_tbl_h10),'#0C7BDC','MarkerFaceAlpha',0.2)
bubblechart(mean(t_tbl_AMM),mean(result_tbl_AMM),var(t_tbl_AMM),'#4B0092','MarkerFaceAlpha',0.2)
bubblelegend('time standard deviation','Location','southeast')
[t1,I1] = sort(mean(t_tbl_opt),'ascend');
[t2,I2] = sort(mean(t_tbl_h10),'ascend');
[t3,I3] = sort(mean(t_tbl_uni),'ascend');
[t4,I4] = sort(mean(t_tbl_AMM),'ascend');
plot(t1,mean(result_tbl_opt(:,I1)),'Color',[255 194 10]/255,'LineWidth', 3)
plot(t3,mean(result_tbl_uni(:,I3)),'Color',[26 255 26]/255,'LineWidth', 3)
plot(t2,mean(result_tbl_h10(:,I2)),'Color',[12 123 220]/255,'LineWidth', 3)
plot(t4,mean(result_tbl_AMM(:,I4)),'Color',[75 0 146]/255,'LineWidth', 3)
xline(std_time,'--','Standard Matrix Multiplication','LineWidth',1.5,'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center');
set(gca, 'YScale', 'log');
%xlim([0,4]);
xlabel('time (s)','FontSize',30);
ylabel('relative error','FontSize',30);
title('Increasing sampling number c: reconstruction performance vs. computation time','fontweight','bold','FontSize',36,'interpreter','latex')
legend('Block AMM','Uniform Sampling','Hutch h=10','Standard AMM','FontSize',20,'interpreter','latex','Location','northeast')
hold off

%Visualization: rolling mean
window_len = 10;
hold on
for i = 1:1:repeat_time
    scatter(t_tbl_opt(i,:),result_tbl_opt(i,:),80,'filled','MarkerFaceColor',[255 194 10]/255,'MarkerEdgeColor',[255 194 10]/255,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.6)
    scatter(t_tbl_uni(i,:),result_tbl_uni(i,:),80,'filled','MarkerFaceColor',[26 255 26]/255,'MarkerEdgeColor',[26 255 26]/255,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.6)
    scatter(t_tbl_h10(i,:),result_tbl_h10(i,:),80,'filled','MarkerFaceColor',[12 123 220]/255,'MarkerEdgeColor',[12 123 220]/255,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.6)
    %scatter(t_tbl_AMM(i,:),result_tbl_AMM(i,:),80,'filled','MarkerFaceColor',[75 0 146]/255,'MarkerEdgeColor',[75 0 146]/255,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.6)
end
% t_opt = t_tbl_opt(:);result_opt = result_tbl_opt(:);
[t1,I1] = sort(t_tbl_opt(:),'ascend');result_opt = result_tbl_opt(:);
[t2,I2] = sort(t_tbl_uni(:),'ascend');result_uni = result_tbl_uni(:);
[t3,I3] = sort(t_tbl_h10(:),'ascend');result_h10 = result_tbl_h10(:);
%[t4,I4] = sort(t_tbl_AMM(:),'ascend');result_AMM = result_tbl_AMM(:);
% [t_temp,I_temp] = sort(t_opt,'ascend');
% temp = cat(1,t_temp',result_opt(I_temp)');
temp1 = cat(1,t1',result_opt(I1)');M1 = movmean(temp1,window_len,2);    
temp2 = cat(1,t2',result_uni(I2)');M2 = movmean(temp2,window_len,2);
temp3 = cat(1,t3',result_h10(I3)');M3 = movmean(temp3,window_len,2);
%temp4 = cat(1,t4',result_AMM(I4)');M4 = movmean(temp4,window_len,2);
plot(M1(1,:),M1(2,:),'Color',[255 194 10]/255,'LineWidth', 3)
plot(M2(1,:),M2(2,:),'Color',[26 255 26]/255,'LineWidth', 3)
plot(M3(1,:),M3(2,:),'Color',[12 123 220]/255,'LineWidth', 3)
%plot(M4(1,:),M4(2,:),'Color',[75 0 146]/255,'LineWidth', 3)
%xline(std_time,'--','Standard Matrix Multiplication','LineWidth',1.5,'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center');
set(gca, 'YScale', 'log');
%xlim([0,1.5]);
xlabel('time (s)','FontSize',30);
ylabel('relative error','FontSize',30);
title('Yellow Taxi data - Feb 2022','fontweight','bold','FontSize',36,'interpreter','latex')
%legend('Block AMM','Uniform Sampling','Hutch h=10','Standard AMM','FontSize',20,'interpreter','latex','Location','northeast')
legend('Block AMM','Uniform Sampling','Hutch h=10','FontSize',20,'interpreter','latex','Location','northeast')
hold off