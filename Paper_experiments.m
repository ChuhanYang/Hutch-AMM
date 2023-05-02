% New experiments with CR line included

%%
%%%Figure 1

% data 1: generate exponential decreasing column data
m = 100;k = 10000;
mu = exp(linspace(50,0,k));
rng(0) 
A = mu + randn(m,k);
%A = rand(m,k);
B = rand(k,m);T = A*B;
block_sz = 100;

% data 2: sparse
m = 10;k = 100000;
A = full(sprand(m,k,0.001));B = exp(4)*full(sprand(k,m,0.01))+ randn(k,m);
T = A*B;
block_sz = 1000;

% data 3: Drift
m = 100;k = 10000;
mu = exp(linspace(5,0,k));
Sigma = eye(k/2);
rng(0) 
A_decay1 = mvnrnd(mu(1:(k/2)),Sigma,m); % A_decay1 = mu(1:(k/2)) + randn(m,k/2);
A_decay2 = mvnrnd(linspace(0,120,k/2),Sigma,m); %A_decay2 = linspace(0,120,k/2) + randn(m,k/2);
A_decay = cat(2,A_decay1,A_decay2);
A = A_decay;
B = rand(k,m);T = A*B;
block_sz = 100;

% data 4: uniform
m = 100;k = 10000;
rng(0) 
A = rand(m,k);
B = rand(k,m);T = A*B;
block_sz = 100;

% partition based on natural sequential
%block_sz = 100;
num_group = ceil(k/block_sz);
index_p = cell(1,num_group); % the indices that each subgroup contains
W = cell(1,num_group);
%idx = randperm(k); % random permutate natural index
for i = 1:1:num_group
    index_p{i} =  ((i-1)*block_sz+1):(i*block_sz) ; %natural sequential
    W{i} = A(:,index_p{i})*B(index_p{i},:);
    p_opt(i) = norm(A(:,index_p{i}))*norm(B(index_p{i},:));
end
p_opt_raw = p_opt;
p_opt = p_opt/sum(p_opt);


repeat_time = 200;
s_list = 10:10:100;
[result_tbl_h1, result_tbl_h10, result_tbl_h50,result_tbl_opt,result_tbl_uni,result_tbl_cr] = deal(zeros(repeat_time,length(s_list)));
for rep_num = 1:1:repeat_time
    disp('repeat round = ')
    disp(rep_num)
    counter = 1;
    for s = s_list
        disp(counter)
        X_h1 = AMM_coarse_hutch_ver4(A,B,s,1,index_p);
        X_h10 = AMM_coarse_hutch_ver4(A,B,s,5,index_p);
        X_h50 = AMM_coarse_hutch_ver4(A,B,s,50,index_p);
        %optimal case
        X_optimal = AMM_true_tracefun_ver2(A,B,s,index_p);
        %uniform sampling case
        X_uniform = AMM_coarse_uni_ver2(A,B,s,index_p);
        X_cr = AMM_true_CR(A,B,s,index_p);
        result_tbl_h1(rep_num,counter) = norm(X_h1-T,'fro')/norm(T,'fro');
        result_tbl_h10(rep_num,counter) = norm(X_h10-T,'fro')/norm(T,'fro');
        result_tbl_h50(rep_num,counter) = norm(X_h50-T,'fro')/norm(T,'fro');
        result_tbl_opt(rep_num,counter) = norm(X_optimal-T,'fro')/norm(T,'fro');
        result_tbl_uni(rep_num,counter) = norm(X_uniform-T,'fro')/norm(T,'fro');
        result_tbl_cr(rep_num,counter) = norm(X_cr-T,'fro')/norm(T,'fro');
        counter = counter+1;
    end
end

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
plot(s_list,mean(result_tbl_cr),'*-','LineWidth', 5,'MarkerSize',20,'Color','k')
set(gca, 'YScale', 'log','FontSize',36)
legend('Hutch AMM h = 1','Hutch AMM h = 5','Hutch AMM h = 50','Optimal AMM','Uniform Sampling','CR','FontSize',36,'interpreter','latex')
hold off;

%%
%%% Figure 3
% data 1: generate exponential decreasing column data
m = 1000;k = 100000;
mu = exp(linspace(50,0,k));
rng(0) 
A = mu + randn(m,k);
%A = rand(m,k);
B = rand(k,m);T = A*B;
block_sz = 1000;

m = 100;k = 10000;
mu = exp(linspace(50,0,k));
rng(0) 
A = mu + randn(m,k);
%A = rand(m,k);
B = rand(k,m);T = A*B;
block_sz = 100;

% data 2: sparse
m = 10;k = 100000;
A = full(sprand(m,k,0.001));B = exp(4)*full(sprand(k,m,0.01))+ randn(k,m);
T = A*B;
block_sz = 1000;

% data 3: Drift
m = 100;k = 10000;
mu = exp(linspace(5,0,k));
Sigma = eye(k/2);
rng(0) 
A_decay1 = mvnrnd(mu(1:(k/2)),Sigma,m); % A_decay1 = mu(1:(k/2)) + randn(m,k/2);
A_decay2 = mvnrnd(linspace(0,120,k/2),Sigma,m); %A_decay2 = linspace(0,120,k/2) + randn(m,k/2);
A_decay = cat(2,A_decay1,A_decay2);
A = A_decay;
B = rand(k,m);T = A*B;
block_sz = 100;

% data 4: uniform
m = 100;k = 10000;
rng(0) 
A = rand(m,k);
B = rand(k,m);T = A*B;
block_sz = 100;


num_group = ceil(k/block_sz);
%p_uniform = ones(1,num_group)/num_group;
index_p = cell(1,num_group);

for i = 1:1:(num_group-1)
    index_p{i} =  ((i-1)*block_sz+1):(i*block_sz) ; %natural sequential
end
index_p{num_group} = ((num_group-1)*block_sz+1):k;

%repeat_time = 20;
repeat_time = 10;


%s_list_uni = floor(linspace(10,400,10)); % TF-IDF
s_list_uni = floor(linspace(10,60000,10)); % yellow taxi
% s_list_uni = floor(linspace(100,2000,10)); % cr sparse
% s_list_uni = floor(linspace(10,2000,10)); % cr exp
% s_list_uni = floor(linspace(10,400,10)); % cr small size synthetic data
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
s_list_opt = floor(linspace(10,30000,10)); % yellow taxi
%s_list_opt = floor(linspace(100,1000,10)); % cr sparse
% s_list_opt = floor(linspace(1000,5000,5));
% s_list_opt = floor(linspace(10,1000,10)); % cr exp
% s_list_opt = floor(linspace(10,300,10)); % cr small size synthetic data
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
s_list_cr = floor(linspace(10,30000,10)); % yellow taxi
%s_list_opt = floor(linspace(1000,400000,10));
%s_list_cr = floor(linspace(100,1000,10)); % cr sparse
%s_list_cr = floor(linspace(1000,5000,5));
% s_list_cr = floor(linspace(10,1000,10)); % cr exp
%s_list_cr = floor(linspace(10,300,10)); % cr small size synthetic data
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

%s_list_h10 = floor(linspace(10,300,10)); % TF-IDF
s_list_h10 = floor(linspace(10,30000,10)); % yellow taxi
%s_list_h10 = floor(linspace(100,1000,10)); % cr sparse
%s_list_h10 = floor(linspace(1000,5000,5));
% s_list_h10 = floor(linspace(10,1000,10)); % cr exp
% s_list_h10 = floor(linspace(10,300,10)); % cr small size synthetic data
result_tbl_h10 = zeros(repeat_time,length(s_list_h10));
t_tbl_h10 = zeros(repeat_time,length(s_list_h10));
for rep_num = 1:1:repeat_time
    disp('h10 repeat round = ')
    disp(rep_num)
    counter = 1;
    for s = s_list_h10 
        disp(counter)
        tic
        X_h10 = AMM_coarse_hutch_ver4(A,B,s,5,index_p);
        t_tbl_h10(rep_num,counter)= toc;
        result_tbl_h10(rep_num,counter) = norm(X_h10-T,'fro')/norm(T,'fro');
        counter = counter+1;
    end
end

% Rolling mean - formal - synthetic
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
xlim([0,1.5]);
xlabel('\textbf{time (s)}','FontSize',46,'interpreter','latex');
ylabel('\textbf{relative error}','FontSize',46,'interpreter','latex');
%title('\textbf{Yellow Taxi Trip, Feb 2022}','FontSize',52,'interpreter','latex')
%legend('Block AMM','Uniform Sampling','Hutch h=10','Standard AMM','FontSize',20,'interpreter','latex','Location','northeast')
legend('Optimal AMM','Uniform Sampling','Hutch AMM h=5','CR','FontSize',36,'interpreter','latex','Location','northeast')
hold off