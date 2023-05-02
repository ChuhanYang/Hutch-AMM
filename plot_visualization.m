% This code is for generating and tuning visualization only

% Visualization: Estimation error comparison

f = figure();
result_tbl = [mean(result_tbl_h1); mean(result_tbl_h10);mean(result_tbl_h50)];
plot(s_list,result_tbl,'o-','LineWidth', 5,'MarkerSize',20)
%ylim([0,0.2])
xlabel('\textbf{AMM sample number}','FontSize',50,'interpreter','latex');
ylabel('\textbf{mean relative error}','FontSize',50,'interpreter','latex');
%title('\textbf{Uniform, Estimation error comparison}','FontSize',80,'interpreter','latex')
hold on;
plot(s_list,mean(result_tbl_opt),'*-','LineWidth', 5,'MarkerSize',20,'Color','#4B0092')
plot(s_list,mean(result_tbl_uni),'*-','LineWidth', 5,'MarkerSize',20,'Color','#40B0A6')
set(gca, 'YScale', 'log','FontSize',36)
legend('Hutch AMM h = 1','Hutch AMM h = 5','Hutch AMM h = 50','Optimal AMM','Uniform Sampling','FontSize',36,'interpreter','latex')
hold off;

%%
%Plot 10 h line with same color but different shade, ver 1
% % Save the result_tbl variable
% save('result_tbl.mat', 'result_tbl');
% f = figure();
% hold on;
% 
% s_list = 10:10:100;
% h_values = 1:10;
% % colors = autumn(length(h_values));
% 
% % Define number of lines and colors
% num_lines = length(h_values);
% dark_color = [0, 0.4470, 0.7410]; % Starting color
% light_color = [0.8500, 0.3250, 0.0980]; % Ending color
% non_diff_color = [0 0 0];
% 
% % Generate custom colormap
% custom_map = zeros(num_lines, 3);
% % for i = 1:num_lines
% %     custom_map(i,:) = dark_color + (light_color - dark_color) * (i-1) / (num_lines-1);
% % end
% % Create colormap with large color steps for i <= 4
% for i = 1:num_lines
%     if i <= 4
%         custom_map(i,:) = dark_color + (light_color - dark_color) * (i-1) / 3;
%     else
%         custom_map(i,:) = non_diff_color;
%     end
% end
% custom_map(10,:) = light_color;
% colormap(custom_map);
% custom_map = [custom_map flipud(linspace(0.9, 0.1, 10))'];
% 
% for h_index = 1:length(h_values)
%     mean_result_tbl = squeeze(mean(result_tbl(:,:,h_index)));
%     %plot(s_list, mean_result_tbl, 'LineWidth', 5, 'Color', colors(h_index,:))
%     plot(s_list, mean_result_tbl, 'LineWidth', 5, 'Color', custom_map(h_index,:))
% 
% %     % Add label to the top of each line
% %     if h_index < 4
% %     text(s_list(1), mean_result_tbl(1), ['Hutch AMM h = ', num2str(h_index)], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 24,'interpreter','latex');
% %     end
% % 
% %     if h_index == 4
% %         text(s_list(5), mean_result_tbl(5), ['Hutch AMM h $\geq$ ', num2str(h_index)], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 24,'interpreter','latex');
% %     end
%     
% end
% 
% 
% 
% xlabel('\textbf{AMM sample number}','FontSize',50,'interpreter','latex');
% ylabel('\textbf{mean relative error}','FontSize',50,'interpreter','latex');
% title('\textbf{Exponential Decreasing, Estimation error comparison}','FontSize',80,'interpreter','latex')
% 
% set(gca, 'YScale', 'log','FontSize',36)
% legend('Hutch AMM h = 1','Hutch AMM h = 2','Hutch AMM h = 3','Hutch AMM h = 4','Hutch AMM h >= 5','FontSize',36,'interpreter','latex')
% hold off;
% % Load the result_tbl variable
% load('result_tbl.mat');

%Plot 10 h line with same color but different shade, ver final?
% Save the result_tbl variable
save('result_tbl.mat', 'result_tbl');
f = figure();
hold on;

s_list = 10:10:100;
h_values = 1:10;
% colors = autumn(length(h_values));

% Define number of lines and colors
num_lines = length(h_values);
color1 = [0, 0.4470, 0.7410]; % color for h=1
color2 = [75,0,146]/255; %color2 = [211 95 183]/255;
color3 = [0,0,0]/255;% color3 = [254 254 98]/255;
color4 = [64 176 166]/255;% color4 = [64 176 166]/255;
non_diff_color  = [0.8500, 0.3250, 0.0980]; % color for h=5

% Generate custom colormap
custom_map = zeros(num_lines, 3);
for i = 1:num_lines
    if i == 1
        custom_map(i,:) = color1;
    elseif i == 2
        custom_map(i,:) = color2;
    elseif i == 3
        custom_map(i,:) = color3;
    elseif i == 4
        custom_map(i,:) = color4;
    else
        custom_map(i,:) = non_diff_color;
    end
end
%custom_map(10,:) = light_color;
colormap(custom_map);
%custom_map = [custom_map flipud(linspace(0.9, 0.1, 10))'];
custom_map = [custom_map 0.6*ones(10,1)];

for h_index = 1:length(h_values)
    mean_result_tbl = squeeze(mean(result_tbl(:,:,h_index)));
    %plot(s_list, mean_result_tbl, 'LineWidth', 5, 'Color', colors(h_index,:))
    plot(s_list, mean_result_tbl, 'LineWidth', 5, 'Color', custom_map(h_index,:))
    
end



xlabel('\textbf{number of samples, c}','FontSize',50,'interpreter','latex');
ylabel('\textbf{mean relative error}','FontSize',50,'interpreter','latex');
%title('\textbf{Exponential Decreasing, Estimation error comparison}','FontSize',80,'interpreter','latex')

set(gca, 'YScale', 'log','FontSize',36)
legend('Hutch AMM h = 1','Hutch AMM h = 2','Hutch AMM h = 3','Hutch AMM h = 4','Hutch AMM h $\geq$ 5','FontSize',36,'interpreter','latex','Location', 'east')
hold off;
% Load the result_tbl variable
load('result_tbl.mat');

%%
%%% Plot shaded area for specified h
h_index = 5;
mean_result_tbl = squeeze(mean(result_tbl(:,:,h_index)));
% Compute 95th and 5th quantiles
quant_95 = prctile(result_tbl(:,:,h_index), 95);
quant_5 = prctile(result_tbl(:,:,h_index), 5);
% Plot line and shaded area
plot(s_list, mean_result_tbl, 'LineWidth', 5, 'Color', 'red');
hold on;
fill([s_list, fliplr(s_list)], [quant_95, fliplr(quant_5)], 'blue', 'FaceAlpha', 0.2, 'LineStyle', 'none');
set(gca, 'YScale', 'log','FontSize',36)
legend('Hutch AMM h = 5', '95% - 5% quantile', 'FontSize', 14);
% Load the result_tbl variable
load('result_tbl.mat');

%%% Plot shaded area for 2 specified h
% Compute 5th and 95th percentiles for h=5 and h=10
h_a = 1;
h_b = 5;
pct5_h5 = squeeze(quantile(result_tbl(:,:,h_a), 0.05));
pct95_h5 = squeeze(quantile(result_tbl(:,:,h_a), 0.95));
pct5_h10 = squeeze(quantile(result_tbl(:,:,h_b), 0.05));
pct95_h10 = squeeze(quantile(result_tbl(:,:,h_b), 0.95));
% Plot mean values for h=5 and h=10 as separate lines
hold on;
plot(s_list, squeeze(mean(result_tbl(:,:,h_a))), 'LineWidth', 5, 'Color', color1)
plot(s_list, squeeze(mean(result_tbl(:,:,h_b))), 'LineWidth', 5, 'Color', non_diff_color)
% Fill in shaded areas with percentiles
x = [s_list, fliplr(s_list)];
y_h5 = [pct5_h5, fliplr(pct95_h5)];
y_h10 = [pct5_h10, fliplr(pct95_h10)];
fill(x, y_h5, 'b', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
fill(x, y_h10, 'r', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
xlabel('\textbf{AMM sample number}','FontSize',50,'interpreter','latex');
ylabel('\textbf{mean relative error}','FontSize',50,'interpreter','latex');
%title('\textbf{Exponential Decreasing, Estimation error comparison}','FontSize',80,'interpreter','latex')
set(gca, 'YScale', 'log','FontSize',36)
%ylim([0,0.7])
legend('Hutch AMM h = 1','Hutch AMM h = 5','FontSize',36,'interpreter','latex','Location', 'best')
hold off;
% Load the result_tbl variable
load('result_tbl.mat');

% mean, add median line
h_a = 1;
h_b = 5;
pct5_h5 = squeeze(quantile(result_tbl(:,:,h_a), 0.05));
pct95_h5 = squeeze(quantile(result_tbl(:,:,h_a), 0.95));
pct5_h10 = squeeze(quantile(result_tbl(:,:,h_b), 0.05));
pct95_h10 = squeeze(quantile(result_tbl(:,:,h_b), 0.95));
% Plot mean values for h=5 and h=10 as separate lines
hold on;
load('result_tbl.mat');
plot(s_list, squeeze(mean(result_tbl(:,:,h_a))), 'LineWidth', 5, 'Color', color1)
plot(s_list, squeeze(mean(result_tbl(:,:,h_b))), 'LineWidth', 5, 'Color', non_diff_color)
% Add dashed lines for median points
load('result_tbl.mat');
y_med_h5 = squeeze(median(result_tbl(:,:,h_a)));
y_med_h10 = squeeze(median(result_tbl(:,:,h_b)));
plot(s_list, y_med_h5, '*-', 'LineWidth', 3, 'Color', color1)
plot(s_list, y_med_h10, '*-', 'LineWidth', 3, 'Color', non_diff_color)
% Fill in shaded areas with percentiles
x = [s_list, fliplr(s_list)];
y_h5 = [pct5_h5, fliplr(pct95_h5)];
y_h10 = [pct5_h10, fliplr(pct95_h10)];
fill(x, y_h5, 'b', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
fill(x, y_h10, 'r', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
xlabel('\textbf{AMM sample number}','FontSize',50,'interpreter','latex');
ylabel('\textbf{mean relative error}','FontSize',50,'interpreter','latex');
%title('\textbf{Exponential Decreasing, Estimation error comparison}','FontSize',80,'interpreter','latex')
set(gca, 'YScale', 'log','FontSize',36)
%ylim([0,0.7])
legend('Hutch AMM h = 1','Hutch AMM h = 5','Median h = 1', 'Median h = 5', 'FontSize',36,'interpreter','latex','Location', 'best')
hold off;
% Load the result_tbl variable
load('result_tbl.mat');

% confidence interval ver 1
% Compute 95% confidence intervals for h=5 and h=10
h_a = 1;
h_b = 5;
alpha = 0.05; % set confidence level
sample_size = size(result_tbl,1);
t_crit = tinv(1 - alpha/2, sample_size); % calculate critical t-value
mean_h5 = squeeze(mean(result_tbl(:,:,h_a)));
mean_h10 = squeeze(mean(result_tbl(:,:,h_b)));
stderr_h5 = squeeze(std(result_tbl(:,:,h_a))) / sqrt(sample_size);
stderr_h10 = squeeze(std(result_tbl(:,:,h_b))) / sqrt(sample_size);
ci_lower_h5 = mean_h5 - t_crit * stderr_h5;
ci_upper_h5 = mean_h5 + t_crit * stderr_h5;
ci_lower_h10 = mean_h10 - t_crit * stderr_h10;
ci_upper_h10 = mean_h10 + t_crit * stderr_h10;

% Plot mean values for h=5 and h=10 as separate lines
hold on;
plot(s_list, mean_h5, 'LineWidth', 5, 'Color', color1)
plot(s_list, mean_h10, 'LineWidth', 5, 'Color', non_diff_color)

% Fill in shaded areas with confidence intervals
x = [s_list, fliplr(s_list)];
y_h5 = [ci_lower_h5, fliplr(ci_upper_h5)];
y_h10 = [ci_lower_h10, fliplr(ci_upper_h10)];
fill(x, y_h5, 'b', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
fill(x, y_h10, 'r', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);

xlabel('\textbf{number of samples, c}','FontSize',50,'interpreter','latex');
ylabel('\textbf{mean relative error}','FontSize',50,'interpreter','latex');
set(gca, 'YScale', 'log','FontSize',36)
legend('Hutch AMM h = 1','Hutch AMM h = 5','FontSize',36,'interpreter','latex','Location', 'best')
hold off;

% Load the result_tbl variable
load('result_tbl.mat');


% confidence interval ver 2 (add median)
h_a = 1;
h_b = 5;
load('result_tbl.mat');
% Calculate mean and standard deviation for each value of s_list
mu_h5 = squeeze(mean(result_tbl(:,:,h_a)));
sigma_h5 = squeeze(std(result_tbl(:,:,h_a)));
mu_h10 = squeeze(mean(result_tbl(:,:,h_b)));
sigma_h10 = squeeze(std(result_tbl(:,:,h_b)));
% Calculate 95% confidence intervals
alpha = 0.05;
z = norminv(1-alpha/2);
ci_h5 = z * sigma_h5 / sqrt(size(result_tbl,2));
ci_h10 = z * sigma_h10 / sqrt(size(result_tbl,2));
% Plot mean values for h=5 and h=10 as separate lines
hold on;
plot(s_list, mu_h5, 'LineWidth', 5, 'Color', color1)
plot(s_list, mu_h10, 'LineWidth', 5, 'Color', non_diff_color)
% Add dashed lines for median points
y_med_h5 = squeeze(median(result_tbl(:,:,h_a)));
y_med_h10 = squeeze(median(result_tbl(:,:,h_b)));
plot(s_list, y_med_h5, '*-', 'LineWidth', 3, 'Color', color1)
plot(s_list, y_med_h10, '*-', 'LineWidth', 3, 'Color', non_diff_color)
% Fill in shaded areas with confidence intervals
x = [s_list, fliplr(s_list)];
y_h5 = [mu_h5-ci_h5, fliplr(mu_h5+ci_h5)];
y_h10 = [mu_h10-ci_h10, fliplr(mu_h10+ci_h10)];
fill(x, y_h5, 'b', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
fill(x, y_h10, 'r', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
xlabel('\textbf{AMM sample number}','FontSize',50,'interpreter','latex');
ylabel('\textbf{mean relative error}','FontSize',50,'interpreter','latex');
set(gca, 'YScale', 'log','FontSize',36)
legend('Hutch AMM h = 1','Hutch AMM h = 5','Median h = 1', 'Median h = 5', 'FontSize',36,'interpreter','latex','Location', 'best')
hold off;

% plus and minus 1 standard dev
h_a = 1;
h_b = 5;
% alpha = 0.05; % set confidence level
sample_size = size(result_tbl,1);
% t_crit = tinv(1 - alpha/2, sample_size); % calculate critical t-value
mean_h5 = squeeze(mean(result_tbl(:,:,h_a)));
mean_h10 = squeeze(mean(result_tbl(:,:,h_b)));
stderr_h5 = squeeze(std(result_tbl(:,:,h_a))) / sqrt(sample_size);
stderr_h10 = squeeze(std(result_tbl(:,:,h_b))) / sqrt(sample_size);
ci_lower_h5 = mean_h5 -  stderr_h5;
ci_upper_h5 = mean_h5 +  stderr_h5;
ci_lower_h10 = mean_h10 - stderr_h10;
ci_upper_h10 = mean_h10 +  stderr_h10;

% Plot mean values for h=5 and h=10 as separate lines
hold on;
plot(s_list, mean_h5, 'LineWidth', 5, 'Color', color1)
plot(s_list, mean_h10, 'LineWidth', 5, 'Color', non_diff_color)

% Fill in shaded areas with confidence intervals
x = [s_list, fliplr(s_list)];
y_h5 = [ci_lower_h5, fliplr(ci_upper_h5)];
y_h10 = [ci_lower_h10, fliplr(ci_upper_h10)];
fill(x, y_h5, 'b', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
fill(x, y_h10, 'r', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);

xlabel('\textbf{number of samples, c}','FontSize',50,'interpreter','latex');
ylabel('\textbf{mean relative error}','FontSize',50,'interpreter','latex');
set(gca, 'YScale', 'log','FontSize',36)
legend('Hutch AMM h = 1','Hutch AMM h = 5','FontSize',36,'interpreter','latex')
hold off;

% Load the result_tbl variable
load('result_tbl.mat');

%%
% Requires R2020a or later
%exportgraphics(f,'Error1.pdf','Resolution',600) 

% Visualization: Estimation error comparison VS time measurements


% f = figure();
% hold on
% for i = 1:1:repeat_time
%     scatter(time_tbl1,result_tbl_opt(i,:),80,'filled','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
%     scatter(time_tbl3,result_tbl_uni(i,:),80,'filled','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
%     scatter(time_tbl2,result_tbl_h10(i,:),80,'filled','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
% end
% [t1,I1] = sort(time_tbl1,'ascend');
% [t2,I2] = sort(time_tbl2,'ascend');
% [t3,I3] = sort(time_tbl3,'ascend');
% plot(t1,mean(result_tbl_opt(:,I1)),'Color','b','LineWidth', 1)
% plot(t3,mean(result_tbl_uni(:,I3)),'Color','r','LineWidth', 1)
% plot(t2,mean(result_tbl_h10(:,I2)),'Color','g','LineWidth', 1)
% set(gca, 'YScale', 'log')
% xlabel('time (s)','FontSize',30);
% ylabel('relative error','FontSize',30);
% title('Increasing sampling number c: reconstruction performance vs. computation time','fontweight','bold','FontSize',36,'interpreter','latex')
% legend('Block AMM','Uniform Sampling','Hutch h=10','FontSize',20,'interpreter','latex','Location','northeast')
% hold off
% 
% % Visualization: Estimation error comparison VS time measurements (3rd version)
% load('workspace_saved/error_vs_time_s500_bsz100_500x1m')
% t_tbl_AMM = t_tbl_AMM(1:5,:);
% result_tbl_AMM = result_tbl_AMM(1:5,:);
% 
% f = figure();
% hold on
% for i = 1:1:repeat_time
%     scatter(t_tbl_opt(i,:),result_tbl_opt(i,:),80,'filled','MarkerFaceColor',[255 194 10]/255,'MarkerEdgeColor',[255 194 10]/255,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.5)
%     scatter(t_tbl_uni(i,:),result_tbl_uni(i,:),80,'filled','MarkerFaceColor',[26 255 26]/255,'MarkerEdgeColor',[26 255 26]/255,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.5)
%     scatter(t_tbl_h10(i,:),result_tbl_h10(i,:),80,'filled','MarkerFaceColor',[12 123 220]/255,'MarkerEdgeColor',[12 123 220]/255,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.5)
%     scatter(t_tbl_AMM(i,:),result_tbl_AMM(i,:),80,'filled','MarkerFaceColor',[75 0 146]/255,'MarkerEdgeColor',[75 0 146]/255,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.5)
% end
% [t1,I1] = sort(mean(t_tbl_opt),'ascend');
% [t2,I2] = sort(mean(t_tbl_h10),'ascend');
% [t3,I3] = sort(mean(t_tbl_uni),'ascend');
% [t4,I4] = sort(mean(t_tbl_AMM),'ascend');
% plot(t1,mean(result_tbl_opt(:,I1)),'Color',[255 194 10]/255,'LineWidth', 3)
% plot(t3,mean(result_tbl_uni(:,I3)),'Color',[26 255 26]/255,'LineWidth', 3)
% plot(t2,mean(result_tbl_h10(:,I2)),'Color',[12 123 220]/255,'LineWidth', 3)
% plot(t4,mean(result_tbl_AMM(:,I4)),'Color',[75 0 146]/255,'LineWidth', 3)
% xline(std_time,'--','Standard Matrix Multiplication','LineWidth',1.5,'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center'); %11.3601202334
% set(gca, 'YScale', 'log')
% %set(gca, 'XScale', 'log')
% xlim([0,4]);
% ylim([0,1.1]);
% xlabel('time (s)','FontSize',30);
% ylabel('relative error','FontSize',30);
% title('Increasing sampling number c: reconstruction performance vs. computation time','fontweight','bold','FontSize',36,'interpreter','latex')
% legend('Block AMM','Uniform Sampling','Hutch h=10','Standard AMM','FontSize',20,'interpreter','latex','Location','northeast')
% hold off
% 
% %Visualization: Bubble chart
% t = tiledlayout(1,1);
% nexttile
% bubblechart(mean(t_tbl_opt),mean(result_tbl_opt),std(t_tbl_opt),'#FFC20A','MarkerFaceAlpha',0.2)
% hold on
% bubblechart(mean(t_tbl_uni),mean(result_tbl_uni),std(t_tbl_uni),'#1AFF1A','MarkerFaceAlpha',0.2)
% bubblechart(mean(t_tbl_h10),mean(result_tbl_h10),std(t_tbl_h10),'#0C7BDC','MarkerFaceAlpha',0.2)
% bubblechart(mean(t_tbl_AMM),mean(result_tbl_AMM),std(t_tbl_AMM),'#4B0092','MarkerFaceAlpha',0.2)
% bubblelegend('time standard deviation','FontSize',14,'Location','southeast');
% % [t1,I1] = sort(mean(t_tbl_opt),'ascend');
% % [t2,I2] = sort(mean(t_tbl_h10),'ascend');
% % [t3,I3] = sort(mean(t_tbl_uni),'ascend');
% % [t4,I4] = sort(mean(t_tbl_AMM),'ascend');
% % plot(t1,mean(result_tbl_opt(:,I1)),'Color',[255 194 10]/255,'LineWidth', 3)
% % plot(t3,mean(result_tbl_uni(:,I3)),'Color',[26 255 26]/255,'LineWidth', 3)
% % plot(t2,mean(result_tbl_h10(:,I2)),'Color',[12 123 220]/255,'LineWidth', 3)
% % plot(t4,mean(result_tbl_AMM(:,I4)),'Color',[75 0 146]/255,'LineWidth', 3)
% xline(std_time,'--','Standard Matrix Multiplication','LineWidth',1.5,'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center');
% set(gca, 'YScale', 'log')
% xlabel('time (s)','FontSize',30);
% ylabel('relative error','FontSize',30);
% title('Bubble chart','fontweight','bold','FontSize',36,'interpreter','latex')
% legend('Block AMM','Uniform Sampling','Hutch h=10','Standard AMM','FontSize',20,'interpreter','latex','Location','northeast')
% hold off
% 
% %Visualization: Bubble chart 2
% t = tiledlayout(1,1);
% nexttile
% hold on
% for i = 1:1:repeat_time
%     scatter(t_tbl_opt(i,:),result_tbl_opt(i,:),80,'filled','+','MarkerFaceColor',[255 194 10]/255,'MarkerEdgeColor',[255 194 10]/255,'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',1)
%     scatter(t_tbl_uni(i,:),result_tbl_uni(i,:),80,'filled','pentagram','MarkerFaceColor',[26 255 26]/255,'MarkerEdgeColor',[26 255 26]/255,'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',1)
%     scatter(t_tbl_h10(i,:),result_tbl_h10(i,:),80,'filled','square','MarkerFaceColor',[12 123 220]/255,'MarkerEdgeColor',[12 123 220]/255,'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',1)
%     scatter(t_tbl_AMM(i,:),result_tbl_AMM(i,:),80,'filled','*','MarkerFaceColor',[75 0 146]/255,'MarkerEdgeColor',[75 0 146]/255,'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',1)
% end
% bubblechart(mean(t_tbl_opt),mean(result_tbl_opt),var(t_tbl_opt),'#FFC20A','MarkerFaceAlpha',0.2)
% bubblechart(mean(t_tbl_uni),mean(result_tbl_uni),var(t_tbl_uni),'#1AFF1A','MarkerFaceAlpha',0.2)
% bubblechart(mean(t_tbl_h10),mean(result_tbl_h10),var(t_tbl_h10),'#0C7BDC','MarkerFaceAlpha',0.2)
% bubblechart(mean(t_tbl_AMM),mean(result_tbl_AMM),var(t_tbl_AMM),'#4B0092','MarkerFaceAlpha',0.2)
% bubblelegend('time standard deviation','Location','southeast')
% [t1,I1] = sort(mean(t_tbl_opt),'ascend');
% [t2,I2] = sort(mean(t_tbl_h10),'ascend');
% [t3,I3] = sort(mean(t_tbl_uni),'ascend');
% [t4,I4] = sort(mean(t_tbl_AMM),'ascend');
% plot(t1,mean(result_tbl_opt(:,I1)),'Color',[255 194 10]/255,'LineWidth', 3)
% plot(t3,mean(result_tbl_uni(:,I3)),'Color',[26 255 26]/255,'LineWidth', 3)
% plot(t2,mean(result_tbl_h10(:,I2)),'Color',[12 123 220]/255,'LineWidth', 3)
% plot(t4,mean(result_tbl_AMM(:,I4)),'Color',[75 0 146]/255,'LineWidth', 3)
% xline(std_time,'--','Standard Matrix Multiplication','LineWidth',1.5,'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center');
% set(gca, 'YScale', 'log');
% %xlim([0,4]);
% xlabel('time (s)','FontSize',30);
% ylabel('relative error','FontSize',30);
% title('Increasing sampling number c: reconstruction performance vs. computation time','fontweight','bold','FontSize',36,'interpreter','latex')
% legend('Block AMM','Uniform Sampling','Hutch h=10','Standard AMM','FontSize',20,'interpreter','latex','Location','northeast')
% hold off
% 
% %Visualization: rolling mean
% window_len = 10;
% hold on
% for i = 1:1:repeat_time
%     scatter(t_tbl_opt(i,:),result_tbl_opt(i,:),80,'filled','MarkerFaceColor',[255 194 10]/255,'MarkerEdgeColor',[255 194 10]/255,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.6)
%     scatter(t_tbl_uni(i,:),result_tbl_uni(i,:),80,'filled','MarkerFaceColor',[26 255 26]/255,'MarkerEdgeColor',[26 255 26]/255,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.6)
%     scatter(t_tbl_h10(i,:),result_tbl_h10(i,:),80,'filled','MarkerFaceColor',[12 123 220]/255,'MarkerEdgeColor',[12 123 220]/255,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.6)
%     %scatter(t_tbl_AMM(i,:),result_tbl_AMM(i,:),80,'filled','MarkerFaceColor',[75 0 146]/255,'MarkerEdgeColor',[75 0 146]/255,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.6)
% end
% % t_opt = t_tbl_opt(:);result_opt = result_tbl_opt(:);
% [t1,I1] = sort(t_tbl_opt(:),'ascend');result_opt = result_tbl_opt(:);
% [t2,I2] = sort(t_tbl_uni(:),'ascend');result_uni = result_tbl_uni(:);
% [t3,I3] = sort(t_tbl_h10(:),'ascend');result_h10 = result_tbl_h10(:);
% %[t4,I4] = sort(t_tbl_AMM(:),'ascend');result_AMM = result_tbl_AMM(:);
% % [t_temp,I_temp] = sort(t_opt,'ascend');
% % temp = cat(1,t_temp',result_opt(I_temp)');
% temp1 = cat(1,t1',result_opt(I1)');M1 = movmean(temp1,window_len,2);    
% temp2 = cat(1,t2',result_uni(I2)');M2 = movmean(temp2,window_len,2);
% temp3 = cat(1,t3',result_h10(I3)');M3 = movmean(temp3,window_len,2);
% %temp4 = cat(1,t4',result_AMM(I4)');M4 = movmean(temp4,window_len,2);
% plot(M1(1,:),M1(2,:),'Color',[255 194 10]/255,'LineWidth', 3)
% plot(M2(1,:),M2(2,:),'Color',[26 255 26]/255,'LineWidth', 3)
% plot(M3(1,:),M3(2,:),'Color',[12 123 220]/255,'LineWidth', 3)
% %plot(M4(1,:),M4(2,:),'Color',[75 0 146]/255,'LineWidth', 3)
% %xline(std_time,'--','Standard Matrix Multiplication','LineWidth',1.5,'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center');
% set(gca, 'YScale', 'log');
% %xlim([0,1.5]);
% xlabel('time (s)','FontSize',30);
% ylabel('relative error','FontSize',30);
% title('Yellow Taxi data - Feb 2022','fontweight','bold','FontSize',36,'interpreter','latex')
% %legend('Block AMM','Uniform Sampling','Hutch h=10','Standard AMM','FontSize',20,'interpreter','latex','Location','northeast')
% legend('Block AMM','Uniform Sampling','Hutch h=10','FontSize',20,'interpreter','latex','Location','northeast')
% hold off


% Rolling mean - formal - synthetic & practical datasets
window_len = 20;
hold on
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
%plot(M1(1,:),M1(2,:),'Color',[255 194 10]/255,'LineWidth', 5)
plot(M1(1,:),M1(2,:),'Color','#4B0092','LineWidth', 5)
%plot(M2(1,:),M2(2,:),'Color',[26 255 26]/255,'LineWidth', 5)
plot(M2(1,:),M2(2,:),'Color','#40B0A6','LineWidth', 5)
plot(M3(1,:),M3(2,:),'Color',[0.8500, 0.3250, 0.0980],'LineWidth', 5)
%plot(M3(1,:),M3(2,:),'Color',[12 123 220]/255,'LineWidth', 5)
%xline(std_time,'--','Standard Matrix Multiplication','LineWidth',1.5,'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center');
for i = 1:1:repeat_time
    %scatter(t_tbl_opt(i,:),result_tbl_opt(i,:),80,'filled','MarkerFaceColor',[255 194 10]/255,'MarkerEdgeColor',[255 194 10]/255,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.6)
    scatter(t_tbl_opt(i,:),result_tbl_opt(i,:),80,'filled','MarkerFaceColor','#4B0092','MarkerEdgeColor','#4B0092','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.6)
    %scatter(t_tbl_uni(i,:),result_tbl_uni(i,:),80,'filled','MarkerFaceColor',[26 255 26]/255,'MarkerEdgeColor',[26 255 26]/255,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.6)
    scatter(t_tbl_uni(i,:),result_tbl_uni(i,:),80,'filled','MarkerFaceColor','#40B0A6','MarkerEdgeColor','#40B0A6','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.6)
    scatter(t_tbl_h10(i,:),result_tbl_h10(i,:),80,'filled','MarkerFaceColor',[0.8500, 0.3250, 0.0980],'MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.6)
    %scatter(t_tbl_h10(i,:),result_tbl_h10(i,:),80,'filled','MarkerFaceColor',[12 123 220]/255,'MarkerEdgeColor',[12 123 220]/255,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.6)
end
set(gca, 'YScale', 'log','FontSize',36);
xlim([0,1.3]);
xlabel('\textbf{time (s)}','FontSize',46,'interpreter','latex');
ylabel('\textbf{relative error}','FontSize',46,'interpreter','latex');
%title('\textbf{Yellow Taxi Trip, Feb 2022}','FontSize',52,'interpreter','latex')
%legend('Block AMM','Uniform Sampling','Hutch h=10','Standard AMM','FontSize',20,'interpreter','latex','Location','northeast')
legend('Optimal AMM','Uniform Sampling','Hutch AMM h=5','FontSize',36,'interpreter','latex','Location','northeast')
hold off