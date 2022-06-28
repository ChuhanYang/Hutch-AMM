% comparison between (Block-) Hutch-AMM & biased p_bar Hutch-AMM & uniform
% sampled AMM

% generate exponential decreasing column data
m = 100;k = 10000;
mu = exp(linspace(5,0,k));
Sigma = eye(k);
rng(0) 
A_decay = mvnrnd(mu,Sigma,m);
A = A_decay;
A = A_decay(:,randperm(k));
B = rand(k,m);T = A*B;

% generate linear increasing
m = 100;k = 10000;
mu = 0:0.001:9.999;
Sigma = eye(k);
rng(0) 
A_decay = mvnrnd(mu,Sigma,m);
A = A_decay;
%A = A_decay(:,randperm(k));
B = 5*rand(k,m);T = A*B;

% generate drift data (exponential decreasing + linear increasing)
m = 100;k = 10000;
mu = exp(linspace(5,0,10000));
Sigma = eye(k/2);
rng(0) 
A_decay1 = mvnrnd(mu(1:(k/2)),Sigma,m);
A_decay2 = mvnrnd(linspace(0,120,k/2),Sigma,m);
A_decay = cat(2,A_decay1,A_decay2);
A = A_decay;
%A = A_decay(:,randperm(k));
B = rand(k,m);T = A*B;

% partition based on natural sequential
block_sz = 100;
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

num_repeat = 100;% number of repeated test
hutch_AMM = zeros(10,num_repeat);
hutch_AMM_bar = zeros(10,num_repeat);
uniform_AMM = zeros(10,num_repeat);
opt_AMM = zeros(1,num_repeat);
v = 3*(1:10);
s = 50; % sample number of AMM
for i = 1:num_repeat 
    disp(i)
    % calculate AMM estimation with optimal probability
    optIndSample_list = randsample(num_group,s,true,p_opt);
    Tr2 = zeros(size(A,1),size(B,2));
    for j = 1:1:s
        t = optIndSample_list(j);
        Tr2 = Tr2 + W{t}/p_opt(t);
    end
    Tr_Opt = Tr2/s;
    opt_AMM(i) = norm(Tr_Opt-T,'fro')/norm(T,'fro');
    %opt_AMM(i) = norm(Tr_Opt-T,'fro')/(norm(A,'fro')*norm(B,'fro'));
    for h = 1:10 % determines sample number of Hutch estimation
        % calculate hutch estimator probability
        p_uniform = ones(1,num_group)/num_group;
        p = zeros(1,num_group);       
        for j = 1:1:num_group
            p(j) = simple_hutchinson(transpose(W{j})*W{j},v(h));
            p(j) = sqrt(p(j));
        end
        p_bar_raw = p;
        p = p/sum(p);
        % run AMM with hutch probability
        hutchIndSample_list = randsample(num_group,s,true,p);
        hutchPlusIndSample_list = randsample(num_group,s,true,p_uniform);
        Tr = zeros(size(A,1),size(B,2));
        Tr_bar = zeros(size(A,1),size(B,2));
        Tr_uniform = zeros(size(A,1),size(B,2));
        sampledset = zeros(1,s);
        for j = 1:1:s
            t = hutchIndSample_list(j);
            t_uniform = hutchPlusIndSample_list(j);            
            t_bar = randsample(num_group,1,true,p_bar_raw/sum(p_bar_raw));                                   
            Tr = Tr + W{t}/p(t);
            Tr_uniform = Tr_uniform + W{t_uniform}/p_uniform(t_uniform);
            % replace the estimated norm with sampled true norm
            if ismember(t_bar,sampledset) ~= true
                p_bar_raw(t_bar) = norm(W{t_bar},'fro');
            end
            sampledset(j) = t_bar;
            Tr_bar = Tr_bar + W{t_bar}/(p_bar_raw(t_bar)/sum(p_bar_raw));                     
        end
        Tr_Hutch = Tr/s;
        Tr_Hutch_bar = Tr_bar/s;
        Tr_Hutch_uniform = Tr_uniform/s; 
        hutch_AMM(h,i) = norm(Tr_Hutch-T,'fro')/norm(T,'fro');
        %hutch_AMM(h,i) = norm(Tr_Hutch-T,'fro')/(norm(A,'fro')*norm(B,'fro'));
        hutch_AMM_bar(h,i) = norm(Tr_Hutch_bar-T,'fro')/norm(T,'fro');
        %hutch_AMM_bar(h,i) = norm(Tr_Hutch_bar-T,'fro')/(norm(A,'fro')*norm(B,'fro'));
        uniform_AMM(h,i) = norm(Tr_Hutch_uniform-T,'fro')/norm(T,'fro');
        %uniform_AMM(h,i) = norm(Tr_Hutch_uniform-T,'fro')/(norm(A,'fro')*norm(B,'fro'));
    end
end

% experiment: fix block zise, increase sample number s, comparing uniform
% sample/Hutch AMM with diff h/Standard AMM
block_sz = 10;
num_group = ceil(k/block_sz);
p_uniform = ones(1,num_group)/num_group;
index_p = cell(1,num_group);
W = cell(1,num_group);
for i = 1:1:num_group
    index_p{i} =  ((i-1)*block_sz+1):(i*block_sz) ; %natural sequential
    W{i} = A(:,index_p{i})*B(index_p{i},:);

end

%s_list = 50:50:500;
s_list = 100:10:490;
result_tbl = zeros(6,length(s_list));
counter = 1;
for s = s_list 
    disp(counter)
    X_h5 = AMM_coarse_hutch(A,B,s,5,index_p);
    X_h10 = AMM_coarse_hutch(A,B,s,10,index_p);
    X_h20 = AMM_coarse_hutch(A,B,s,20,index_p);
    X_h50 = AMM_coarse_hutch(A,B,s,50,index_p);
    %optimal case
    X_optimal = AMM_true(A,B,s,index_p);

    %uniform sampling case
    X_uniform = zeros(size(A,1),size(B,2));
    UniformSample_list = randsample(num_group,s,true,p_uniform);
    for j = 1:1:s
        t_uniform = UniformSample_list(j);           
        X_uniform = X_uniform + W{t_uniform}/p_uniform(t_uniform);                    
    end
    X_uniform = X_uniform/s;
    result_tbl(1,counter) = norm(X_h5-T,'fro')/norm(T,'fro');
    result_tbl(2,counter) = norm(X_h10-T,'fro')/norm(T,'fro');
    result_tbl(3,counter) = norm(X_h20-T,'fro')/norm(T,'fro');
    result_tbl(4,counter) = norm(X_h50-T,'fro')/norm(T,'fro');
    result_tbl(5,counter) = norm(X_optimal-T,'fro')/norm(T,'fro');
    result_tbl(6,counter) = norm(X_uniform-T,'fro')/norm(T,'fro');
    counter = counter+1;
end

plot(s_list,result_tbl,'o-')
%ylim([0,0.18])
xlabel('AMM sample number','FontSize',15);
ylabel('relative error','FontSize',15);
title('Exponential decreasing, Estimation error comparison, block size = 10','FontSize',20)
legend('Hutch AMM h = 5','Hutch AMM h = 10','Hutch AMM h = 20','Hutch AMM h = 50','Optimal AMM','Uniform Sampling')

% visualize mean values:
block_sz = 10;
num_group = ceil(k/block_sz);
p_uniform = ones(1,num_group)/num_group;
index_p = cell(1,num_group);
W = cell(1,num_group);
for i = 1:1:num_group
    index_p{i} =  ((i-1)*block_sz+1):(i*block_sz) ; %natural sequential
    W{i} = A(:,index_p{i})*B(index_p{i},:);

end

repeat_time = 20;
s_list = 100:10:490;
[result_tbl_h5, result_tbl_h10, result_tbl_h20,result_tbl_h50,result_tbl_opt,result_tbl_uni] = deal(zeros(repeat_time,length(s_list)));
for rep_num = 1:1:repeat_time
    disp('repeat round = ')
    disp(rep_num)
    counter = 1;
    for s = s_list 
        disp(counter)
        X_h5 = AMM_coarse_hutch(A,B,s,5,index_p);
        X_h10 = AMM_coarse_hutch(A,B,s,10,index_p);
        X_h20 = AMM_coarse_hutch(A,B,s,20,index_p);
        X_h50 = AMM_coarse_hutch(A,B,s,50,index_p);
        %optimal case
        X_optimal = AMM_true(A,B,s,index_p);

        %uniform sampling case
        X_uniform = zeros(size(A,1),size(B,2));
        UniformSample_list = randsample(num_group,s,true,p_uniform);
        for j = 1:1:s
            t_uniform = UniformSample_list(j);           
            X_uniform = X_uniform + W{t_uniform}/p_uniform(t_uniform);                    
        end
        X_uniform = X_uniform/s;
        result_tbl_h5(rep_num,counter) = norm(X_h5-T,'fro')/norm(T,'fro');
        result_tbl_h10(rep_num,counter) = norm(X_h10-T,'fro')/norm(T,'fro');
        result_tbl_h20(rep_num,counter) = norm(X_h20-T,'fro')/norm(T,'fro');
        result_tbl_h50(rep_num,counter) = norm(X_h50-T,'fro')/norm(T,'fro');
        result_tbl_opt(rep_num,counter) = norm(X_optimal-T,'fro')/norm(T,'fro');
        result_tbl_uni(rep_num,counter) = norm(X_uniform-T,'fro')/norm(T,'fro');
        counter = counter+1;
    end
end

result_tbl = [mean(result_tbl_h5); mean(result_tbl_h10);mean(result_tbl_h20);mean(result_tbl_h50);mean(result_tbl_opt);mean(result_tbl_uni)];
plot(s_list,result_tbl,'o-')
%ylim([0,0.18])
xlabel('AMM sample number','FontSize',15);
ylabel('mean relative error','FontSize',15);
title('Exponential decreasing with permutation, Estimation error comparison, block size = 10','FontSize',20)
legend('Hutch AMM h = 5','Hutch AMM h = 10','Hutch AMM h = 20','Hutch AMM h = 50','Optimal AMM','Uniform Sampling')


% visualization: performance comparison
fig = figure;
row_names = arrayfun(@num2str,v,'uni',0);
row_names = [row_names 'opt'];
subplot(3,1,1);
boxplot([hutch_AMM' opt_AMM'],'Labels',row_names)
%ylim([0.011 0.018])
title('Subplot 1: Hutch-AMM')
subplot(3,1,2);
boxplot([hutch_AMM_bar' opt_AMM'],'Labels',row_names)
%ylim([0.011 0.018])
title('Subplot 2: Biased Hutch-AMM')
subplot(3,1,3);
boxplot([uniform_AMM' opt_AMM'],'Labels',row_names)
%ylim([0.011 0.018])
title('Subplot 3: uniform')
sgtitle('Exponential decreasing with permutation, sample number = 50, block size = 100','FontSize',20)
han=axes(fig,'visible','off');
han.YLabel.Visible='on';
han.XLabel.Visible='on';
ylabel(han,'standard relative error','FontSize',15);
xlabel(han,'number of matrix-vector multiplication','FontSize',15);

% visualization: data property
m = 100;k = 10000;
mu = exp(linspace(5,0,10000));
Sigma = eye(k);
rng(0) 
A_exp = mvnrnd(mu,Sigma,m);
A_permute = A_exp(:,randperm(k));

mu = 0:0.001:9.999;
Sigma = eye(k);
rng(0) 
A_linear = mvnrnd(mu,Sigma,m);

mu = exp(linspace(5,0,10000));
Sigma = eye(k/2);
rng(0) 
A_decay1 = mvnrnd(mu(1:(k/2)),Sigma,m);
A_decay2 = mvnrnd(linspace(0,120,k/2),Sigma,m);
A_drift = cat(2,A_decay1,A_decay2);

mean_v_1 = zeros(1,k);mean_v_2 = zeros(1,k);mean_v_3 = zeros(1,k);mean_v_4 = zeros(1,k);
for i = 1:1:k
    mean_v_1(i) = mean(A_exp(:,i));
    mean_v_2(i) = mean(A_permute(:,i));
    mean_v_3(i) = mean(A_linear(:,i));
    mean_v_4(i) = mean(A_drift(:,i));
end

fig = figure;
subplot(2,2,1);
plot(linspace(1,k,k),mean_v_1)
title('Subplot 1: Exponential decreasing')
subplot(2,2,2);
plot(linspace(1,k,k),mean_v_2)
title('Subplot 2: Exponential decreasing After permutation')
subplot(2,2,3);
plot(linspace(1,k,k),mean_v_3)
title('Subplot 2: Linear increasing')
subplot(2,2,4);
plot(linspace(1,k,k),mean_v_4)
han=axes(fig,'visible','off');
han.YLabel.Visible='on';
han.XLabel.Visible='on';
title('Subplot 2: Drift: exponential decreasing + linear increasing')
ylabel(han,'y-axis: mean value of this column','FontSize',15);
xlabel(han,'x-axis: column index','FontSize',15);
sgtitle('mean value of each column from A')

% visualization: runtime of diff block_sz
N = 5000; % total allowed sampled column number
x = 1:N;
block_sz_List = x(~(rem(N, x))); % find divisor for N as block_sz
time_tbl = zeros(1,20);
counter = 1;
for block_sz = block_sz_List
    disp(counter)
    s = ceil(N/block_sz);
    h = 10;    
    f = @() Hutch_AMM(A, B, s, h, block_sz);
    time_tbl(counter) = timeit(f);
    counter = counter+1;
end

plot(block_sz_List,time_tbl,'o-')
set(gca, 'YScale', 'log')
xlabel('block size');
ylabel('time (s)');
title('Avergae runtime of Bloack Hutch AMM')

% visualization: runtime of AMM and Block-Hutch-AMM
N = 5000; % total allowed sampled column number
x = 1:N;
block_sz_List = x(~(rem(N, x))); % find divisor for N as block_sz
time_tbl1 = zeros(1,20);time_tbl2 = zeros(1,20);
counter = 1;
for block_sz = block_sz_List
    disp(counter)
    num_group = ceil(k/block_sz);
    index_p = cell(1,num_group);
    W = cell(1,num_group);
    for i = 1:1:num_group
        index_p{i} =  ((i-1)*block_sz+1):(i*block_sz) ; %natural sequential
    end
    s = ceil(N/block_sz);
    h = 10;
    f1 = @() AMM_true_tracefun(A,B,s,index_p);
    f2 = @() AMM_coarse_hutch(A,B,s,h,index_p);
    time_tbl1(counter) = timeit(f1);
    time_tbl2(counter) = timeit(f2);
    counter = counter+1;
end


plot(block_sz_List,time_tbl1,'o-',block_sz_List,time_tbl2,'*-')
set(gca, 'YScale', 'log')
xlabel('block size');
ylabel('time (s)');
title('Avergae runtime of Bloack Hutch AMM with predefined groups')
legend('Block AMM','Hutch Block AMM')