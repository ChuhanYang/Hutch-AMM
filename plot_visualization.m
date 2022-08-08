% This code is for generating and tuning visualization only

% Visualization: Estimation error comparison
load('workspace_saved/ExpDec_h10_ErrEst')
load('workspace_saved/LinInc_h10_ErrEst')
load('workspace_saved/Drift_h10_ErrEst')
load('workspace_saved/Uniform_h10_ErrEst')

f = figure();
result_tbl = [mean(result_tbl_h1);mean(result_tbl_h3);mean(result_tbl_h5); mean(result_tbl_h10);mean(result_tbl_h20);mean(result_tbl_h50)];
plot(s_list,result_tbl,'o-','LineWidth', 1,'MarkerSize',10)
%ylim([0,0.18])
xlabel('AMM sample number','FontSize',30);
ylabel('mean relative error','FontSize',30);
title('Block Size = 20, Estimation error comparison','fontweight','bold','FontSize',36,'interpreter','latex')
hold on;
plot(s_list,mean(result_tbl_opt),'*-','LineWidth', 1,'MarkerSize',10)
plot(s_list,mean(result_tbl_uni),'*-','LineWidth', 1,'MarkerSize',10)
legend('Hutch AMM h = 1','Hutch AMM h = 3','Hutch AMM h = 5','Hutch AMM h = 10','Hutch AMM h = 20','Hutch AMM h = 50','Optimal AMM','Uniform Sampling','FontSize',20,'interpreter','latex')
hold off;

% Requires R2020a or later
%exportgraphics(f,'Error1.pdf','Resolution',600) 

% Visualization: Estimation error comparison VS time measurements
load('workspace_saved/error_vs_time')
load('workspace_saved/error_vs_time_bsz100')

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
legend('Block AMM','Hutch h=10','Uniform Sampling','FontSize',20,'interpreter','latex','Location','northeast')
hold off