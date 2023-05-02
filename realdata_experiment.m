% Selected experiments on practical datasets

%% Read data
% unit_change = 1000/3600;
% A = X(:,1:2700)*unit_change;
% rou_jam = 1/7.5;
% c_cong = -15*unit_change;
% c_free = 100*unit_change;
% B = real(rou_jam * (1 + c_free/c_cong * log(1-A/c_free)).^(-1))';
% [m,k] = size(A);T=A*B;
% 
% % Movie property data
% load("X_movie.mat")
% A = double(count_matrix);[m,k] = size(A);B=A';
% T=A*B;
% 
% % Wiki People data
% load("X_wikipeople.mat")
% A = full(wiki_matrix);A = A(randsample(size(A,1),1000,false),:);[m,k] = size(A);B=A';
% T=A*B;
% 
% % 20news data
% load("X_20news.mat")
% % readtable('X_20news.mat')
% % dlmread("X_20news.mat")
% A = full(newsgroup);A = A(randsample(size(A,1),1000,false),:);[m,k] = size(A);B=A';
% T=A*B;
% 
% % recv1 data
% load("recv1.mat")
% A = recv1_matrix;A = full(A(randsample(size(A,1),1000,false),:));[m,k] = size(A);B=A';
% T=A*B;
% 
% % Beijing weather data
% A = readtable('PRSA_Data_Wanshouxigong_20130301-20170228.csv');
% A = PRSADataChangping2013030120170228(:,[6:15,17]);
% A = table2array(A)';
% A = fillmissing(A,'previous',2);
% %A = normr(A);
% [m,k] = size(A);B=A';T=A*B;
% 
% % household electric power consumption data
% A = householdpowerconsumption(:,3:end);
% A = table2array(A)';
% A = fillmissing(A,'previous',2);
% [m,k] = size(A);B=A';T=A*B;


% TF-IDF TDT2 data
load("TDT2.mat")
fea_sample = fea(randsample(size(fea,1),1000,false),:);
fea = tfidf(fea_sample);
A = full(fea);[m,k] = size(A);B = A';
T=A*B;

% yellow taxi data
%A = yellowtaxinum(:,2:end);
A = dfnumonehot(:,2:end);
A = table2array(A)';
[m,k] = size(A);B=A';T=A*B;

%%
block_sz = 100;
num_group = ceil(k/block_sz);
index_p = cell(1,num_group);

for i = 1:1:(num_group-1)
    index_p{i} =  ((i-1)*block_sz+1):(i*block_sz) ; %natural sequential
end
index_p{num_group} = ((num_group-1)*block_sz+1):k;

%repeat_time = 20;
repeat_time = 10;

% s_list_h10 = floor(linspace(10,300,10)); % TF-IDF
% s_list_h10 = floor(linspace(10,1000,10)); % weather
% s_list_h10 = floor(linspace(10,1000,5)); % synthetic MA
% s_list_h10 = floor(linspace(10,7000,5));
% s_list_h10 = floor(linspace(10,20000,5)); % yellow taxi (old)
% s_list_h10 = floor(linspace(10,30000,10)); % yellow taxi
%s_list_h10 = floor(linspace(1000,400000,10));
%s_list_h10 = floor(linspace(100,1000,10)); % cr sparse
s_list_h10 = floor(linspace(1000,5000,5));
result_tbl_h10 = zeros(repeat_time,length(s_list_h10));
t_tbl_h10 = zeros(repeat_time,length(s_list_h10));
for rep_num = 1:1:repeat_time
    disp('h10 repeat round = ')
    disp(rep_num)
    counter = 1;
    for s = s_list_h10 
        disp(counter)
        tic
        X_h10 = AMM_coarse_hutch_ver4(A,B,s,10,index_p);
        t_tbl_h10(rep_num,counter)= toc;
        result_tbl_h10(rep_num,counter) = norm(X_h10-T,'fro')/norm(T,'fro');
        counter = counter+1;
    end
end


% s_list_uni = floor(linspace(10,600,10)); % TF-IDF
% s_list_uni = floor(linspace(10,1000,10)); % weather
% s_list_uni = floor(linspace(10,1000,5)); % synthetic MA
% s_list_uni = floor(linspace(10,7000,5));
% s_list_uni = floor(linspace(10,60000,10)); % yellow taxi
%s_list_uni = floor(linspace(1000,400000,10));
% s_list_uni = floor(linspace(100,2000,10)); % cr sparse
s_list_uni = floor(linspace(1000,10000,5));
result_tbl_uni = zeros(repeat_time,length(s_list_uni));
t_tbl_uni = zeros(repeat_time,length(s_list_uni));
for rep_num = 1:1:repeat_time
    disp('uni repeat round = ')
    disp(rep_num)
    counter = 1;
    for s = s_list_uni 
        disp(counter)
        tic
        %X_uniform = AMM_coarse_uni_ver2(A,B,s,index_p,false); %without replacement
        X_uniform = AMM_coarse_uni_ver2(A,B,s,index_p); %with replacement
        t_tbl_uni(rep_num,counter)= toc;
        result_tbl_uni(rep_num,counter) = norm(X_uniform-T,'fro')/norm(T,'fro');        
        counter = counter+1;
    end
end

%s_list_opt = floor(linspace(10,300,10)); % TF-IDF
% s_list_opt = floor(linspace(10,1000,10)); % weather
% s_list_opt = floor(linspace(10,1000,5)); % synthetic MA
% s_list_opt = floor(linspace(10,7000,5));
% s_list_opt = floor(linspace(10,30000,10)); % yellow taxi
%s_list_opt = floor(linspace(1000,400000,10));
%s_list_opt = floor(linspace(100,1000,10)); % cr sparse
s_list_opt = floor(linspace(1000,5000,5));
result_tbl_opt = zeros(repeat_time,length(s_list_opt));
t_tbl_opt = zeros(repeat_time,length(s_list_opt));
for rep_num = 1:1:repeat_time
    disp('opt repeat round = ')
    disp(rep_num)
    counter = 1;
    for s = s_list_opt 
        disp(counter)
        tic
        X_optimal = AMM_true_tracefun_ver2(A,B,s,index_p);
        t_tbl_opt(rep_num,counter)= toc;
        result_tbl_opt(rep_num,counter) = norm(X_optimal-T,'fro')/norm(T,'fro');
        counter = counter+1;
    end
end


%%%% add CR method
%s_list_cr = floor(linspace(10,300,10)); % TF-IDF
% s_list_opt = floor(linspace(10,1000,10)); % weather
% s_list_opt = floor(linspace(10,1000,5)); % synthetic MA
% s_list_opt = floor(linspace(10,7000,5));
% s_list_cr = floor(linspace(10,30000,10)); % yellow taxi
%s_list_opt = floor(linspace(1000,400000,10));
%s_list_cr = floor(linspace(100,1000,10)); % cr sparse
s_list_cr = floor(linspace(1000,5000,5));
result_tbl_cr = zeros(repeat_time,length(s_list_cr));
t_tbl_cr = zeros(repeat_time,length(s_list_cr));
for rep_num = 1:1:repeat_time
    disp('cr repeat round = ')
    disp(rep_num)
    counter = 1;
    for s = s_list_cr 
        disp(counter)
        tic
        X_cr = AMM_true_CR(A,B,s,index_p);
        t_tbl_cr(rep_num,counter)= toc;
        result_tbl_cr(rep_num,counter) = norm(X_cr-T,'fro')/norm(T,'fro');
        counter = counter+1;
    end
end


%% Rolling mean - visualization
window_len = 20;
hold on
% t_opt = t_tbl_opt(:);result_opt = result_tbl_opt(:);
[t1,I1] = sort(t_tbl_opt(:),'ascend');result_opt = result_tbl_opt(:);
[t2,I2] = sort(t_tbl_uni(:),'ascend');result_uni = result_tbl_uni(:);
[t3,I3] = sort(t_tbl_h10(:),'ascend');result_h10 = result_tbl_h10(:);
[t4,I4] = sort(t_tbl_cr(:),'ascend');result_cr = result_tbl_cr(:);
%[t4,I4] = sort(t_tbl_AMM(:),'ascend');result_AMM = result_tbl_AMM(:);
% [t_temp,I_temp] = sort(t_opt,'ascend');
% temp = cat(1,t_temp',result_opt(I_temp)');
temp1 = cat(1,t1',result_opt(I1)');M1 = movmean(temp1,window_len,2);    
temp2 = cat(1,t2',result_uni(I2)');M2 = movmean(temp2,window_len,2);
temp3 = cat(1,t3',result_h10(I3)');M3 = movmean(temp3,window_len,2);
temp4 = cat(1,t4',result_cr(I4)');M4 = movmean(temp4,window_len,2);
%temp4 = cat(1,t4',result_AMM(I4)');M4 = movmean(temp4,window_len,2);
%plot(M1(1,:),M1(2,:),'Color',[255 194 10]/255,'LineWidth', 5)
plot(M1(1,:),M1(2,:),'Color','#4B0092','LineWidth', 5)
%plot(M2(1,:),M2(2,:),'Color',[26 255 26]/255,'LineWidth', 5)
plot(M2(1,:),M2(2,:),'Color','#40B0A6','LineWidth', 5)
plot(M3(1,:),M3(2,:),'Color',[0.8500, 0.3250, 0.0980],'LineWidth', 5)
%plot(M3(1,:),M3(2,:),'Color',[12 123 220]/255,'LineWidth', 5)
plot(M4(1,:),M4(2,:),'Color','k','LineWidth', 5)
%xline(std_time,'--','Standard Matrix Multiplication','LineWidth',1.5,'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center');
for i = 1:1:repeat_time
    %scatter(t_tbl_opt(i,:),result_tbl_opt(i,:),80,'filled','MarkerFaceColor',[255 194 10]/255,'MarkerEdgeColor',[255 194 10]/255,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.6)
    scatter(t_tbl_opt(i,:),result_tbl_opt(i,:),80,'filled','MarkerFaceColor','#4B0092','MarkerEdgeColor','#4B0092','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.6)
    %scatter(t_tbl_uni(i,:),result_tbl_uni(i,:),80,'filled','MarkerFaceColor',[26 255 26]/255,'MarkerEdgeColor',[26 255 26]/255,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.6)
    scatter(t_tbl_uni(i,:),result_tbl_uni(i,:),80,'filled','MarkerFaceColor','#40B0A6','MarkerEdgeColor','#40B0A6','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.6)
    scatter(t_tbl_h10(i,:),result_tbl_h10(i,:),80,'filled','MarkerFaceColor',[0.8500, 0.3250, 0.0980],'MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.6)
    scatter(t_tbl_cr(i,:),result_tbl_cr(i,:),80,'filled','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.6)
end
set(gca, 'YScale', 'log','FontSize',36);
%xlim([0,0.13]);
xlabel('\textbf{time (s)}','FontSize',46,'interpreter','latex');
ylabel('\textbf{relative error}','FontSize',46,'interpreter','latex');
%title('\textbf{Yellow Taxi Trip, Feb 2022}','FontSize',52,'interpreter','latex')
%legend('Block AMM','Uniform Sampling','Hutch h=10','Standard AMM','FontSize',20,'interpreter','latex','Location','northeast')
legend('Optimal AMM','Uniform Sampling','Hutch AMM h=5','CR','FontSize',36,'interpreter','latex','Location','northeast')
hold off