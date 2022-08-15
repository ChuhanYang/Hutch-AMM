% comparison between (Block-) Hutch-AMM & biased p_bar Hutch-AMM & uniform
% sampled AMM

% generate exponential decreasing column data
%m = 800;k = 10000;
m = 500;k = 1000000;
mu = exp(linspace(50,0,k));
%mu = exp(linspace(5,0,k));
%Sigma = eye(k);
rng(0) 
%A_decay = mvnrnd(mu,Sigma,m);
A_decay = mu + randn(m,k);
A = A_decay;
%A = A_decay(:,randperm(k));
%A = rand(m,k);
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
%block_sz = 5;
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
%%repeat_time = 50;
s_list = 10:10:490;
[result_tbl_h1,result_tbl_h3,result_tbl_h5, result_tbl_h10, result_tbl_h20,result_tbl_h50,result_tbl_opt,result_tbl_uni] = deal(zeros(repeat_time,length(s_list)));
for rep_num = 1:1:repeat_time
    disp('repeat round = ')
    disp(rep_num)
    counter = 1;
    for s = s_list 
        disp(counter)
        X_h1 = AMM_coarse_hutch(A,B,s,1,index_p);
        X_h3 = AMM_coarse_hutch(A,B,s,3,index_p);
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
        result_tbl_h1(rep_num,counter) = norm(X_h1-T,'fro')/norm(T,'fro');
        result_tbl_h3(rep_num,counter) = norm(X_h3-T,'fro')/norm(T,'fro');
        result_tbl_h5(rep_num,counter) = norm(X_h5-T,'fro')/norm(T,'fro');
        result_tbl_h10(rep_num,counter) = norm(X_h10-T,'fro')/norm(T,'fro');
        result_tbl_h20(rep_num,counter) = norm(X_h20-T,'fro')/norm(T,'fro');
        result_tbl_h50(rep_num,counter) = norm(X_h50-T,'fro')/norm(T,'fro');
        result_tbl_opt(rep_num,counter) = norm(X_optimal-T,'fro')/norm(T,'fro');
        result_tbl_uni(rep_num,counter) = norm(X_uniform-T,'fro')/norm(T,'fro');
        counter = counter+1;
    end
end

%save('ExpDec_h10_ErrEst')
%save('LinInc_h10_ErrEst')
%save('Drift_h10_ErrEst')
%save('Uniform_h10_ErrEst')
save('ExpDec_bsz20_ErrEst')

result_tbl = [mean(result_tbl_h1);mean(result_tbl_h3);mean(result_tbl_h5); mean(result_tbl_h10);mean(result_tbl_h20);mean(result_tbl_h50)];
plot(s_list,result_tbl,'o-')
%ylim([0,0.18])
xlabel('AMM sample number','FontSize',15);
ylabel('mean relative error','FontSize',15);
title('Exponential decreasing, Estimation error comparison','FontSize',20)
hold on;
plot(s_list,mean(result_tbl_opt),'*-')
plot(s_list,mean(result_tbl_uni),'*-')
legend('Hutch AMM h = 1','Hutch AMM h = 3','Hutch AMM h = 5','Hutch AMM h = 10','Hutch AMM h = 20','Hutch AMM h = 50','Optimal AMM','Uniform Sampling')
hold off;

% visualization: shaded line plot
fig = gcf;
options.handle = fig;
options.alpha      = 0.5;
options.line_width = 2;
options.error      = 'c95';
options.color_area = [243 169 114]./255;    % Orange theme
options.color_line = [236 112  22]./255;
plot_areaerrorbar(result_tbl_h5,options)
set(gca,'XTick',linspace(1,length(s_list),9))
set(gca,'XTickLabel',linspace(10,490,9))
%ylim([0,1])
hold on
plot_areaerrorbar(result_tbl_h50)
hold on
options.color_area = [0.4 0.4 0.4];    
options.color_line = [0.85 0.85 0.85];% Grey theme
plot_areaerrorbar(result_tbl_uni,options);
hold on
options.color_area = [0 204 102]./255;    
options.color_line = [0 153 75]./255; % Green theme
plot_areaerrorbar(result_tbl_opt,options);
hold on
xlabel('AMM sample number','FontSize',15);
ylabel('mean relative error','FontSize',15);
title('Exponential decreasing with block size = 10,#RepeatTest=50','FontSize',20)
legend('', 'Hutch AMM h = 5', '','Hutch AMM h = 50','','Uniform Sampling','','Optimal AMMA')
hold off


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
N = 2500; % total allowed sampled column number
x = 1:N;
block_sz_List = x(~(rem(N, x))); % find divisor for N as block_sz
%block_sz_List = block_sz_List(5:end); % remove block_sz < 10
%time_tbl1 = zeros(1,20);time_tbl2 = zeros(1,20);time_tbl3 = zeros(1,20);
[time_tbl1,time_tbl2,time_tbl3,time_tbl4,time_tbl5] = deal(zeros(1,length(block_sz_List)));
counter = 1;
for block_sz = block_sz_List
    disp(counter)
    num_group = ceil(k/block_sz);
    index_p = cell(1,num_group);
    for i = 1:1:num_group
        index_p{i} =  ((i-1)*block_sz+1):(i*block_sz) ; %natural sequential
    end
    s = ceil(N/block_sz);
    f1 = @() AMM_true_tracefun_ver2(A,B,s,index_p);
    f2 = @() AMM_coarse_hutch_ver3(A,B,s,10,index_p);
    f3 = @() AMM_coarse_uni_ver2(A,B,s,index_p);
    f4 = @() AMM_coarse_hutch_ver3(A,B,s,1,index_p);
    f5 = @() AMM_coarse_hutch_ver3(A,B,s,50,index_p);
    time_tbl1(counter) = timeit(f1);
    time_tbl2(counter) = timeit(f2);
    time_tbl3(counter) = timeit(f3);
    time_tbl4(counter) = timeit(f4);
    time_tbl5(counter) = timeit(f5);
    counter = counter+1;
end

plot(block_sz_List,time_tbl1,'o-',block_sz_List,time_tbl4,'*-',block_sz_List,time_tbl2,'*-',block_sz_List,time_tbl5,'*-',block_sz_List,time_tbl3,'+-')
set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')
xlabel('block size');
ylabel('time (s)');
title('Size 400 by 10000')
legend('Block AMM','Hutch h=1','Hutch h=10','Hutch h=50','Uniform Sampling')

save('time_400x10k')

% visualization: relative error vs runtime (1st version)
block_sz = 25;
num_group = ceil(k/block_sz);
p_uniform = ones(1,num_group)/num_group;
index_p = cell(1,num_group);
W = cell(1,num_group);
for i = 1:1:num_group
    index_p{i} =  ((i-1)*block_sz+1):(i*block_sz) ; %natural sequential
    W{i} = A(:,index_p{i})*B(index_p{i},:);
end

repeat_time = 10;
s_list = linspace(100,1000,10);
s_list2 = linspace(500,5000,10);
[result_tbl_h1, result_tbl_h10,result_tbl_h50,result_tbl_opt,result_tbl_uni] = deal(zeros(repeat_time,length(s_list)));
for rep_num = 1:1:repeat_time
    disp('repeat round = ')
    disp(rep_num)
    counter = 1;
    for s = s_list 
        disp(counter)
        %X_h1 = AMM_coarse_hutch(A,B,s,1,index_p);
        X_h10 = AMM_coarse_hutch(A,B,s,10,index_p);
        %X_h50 = AMM_coarse_hutch(A,B,s,50,index_p);
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
        %result_tbl_h1(rep_num,counter) = norm(X_h1-T,'fro')/norm(T,'fro');
        result_tbl_h10(rep_num,counter) = norm(X_h10-T,'fro')/norm(T,'fro');
        %result_tbl_h50(rep_num,counter) = norm(X_h50-T,'fro')/norm(T,'fro');
        result_tbl_opt(rep_num,counter) = norm(X_optimal-T,'fro')/norm(T,'fro');
        result_tbl_uni(rep_num,counter) = norm(X_uniform-T,'fro')/norm(T,'fro');
        counter = counter+1;
    end
end


[time_tbl1,time_tbl2,time_tbl3,time_tbl4,time_tbl5] = deal(zeros(1,length(s_list)));
counter = 1;
for s = s_list
    disp(counter)
    num_group = ceil(k/block_sz);
    index_p = cell(1,num_group);
    W = cell(1,num_group);
    for i = 1:1:num_group
        index_p{i} =  ((i-1)*block_sz+1):(i*block_sz) ; %natural sequential
    end
    f1 = @() AMM_true_tracefun(A,B,s,index_p);
    f2 = @() AMM_coarse_hutch(A,B,s,10,index_p);
    f3 = @() AMM_coarse_uni(A,B,s,index_p);
    %f4 = @() AMM_coarse_hutch(A,B,s,1,index_p);
    %f5 = @() AMM_coarse_hutch(A,B,s,50,index_p);
    time_tbl1(counter) = timeit(f1);
    time_tbl2(counter) = timeit(f2);
    time_tbl3(counter) = timeit(f3);
    %time_tbl4(counter) = timeit(f4);
    %time_tbl5(counter) = timeit(f5);
    counter = counter+1;
end

plot(time_tbl1,mean(result_tbl_opt),'o-',time_tbl2,mean(result_tbl_h10),'*-',time_tbl3,mean(result_tbl_uni),'+-')
xlabel('time (s)');
ylabel('relative error');
title('Increasing sampling number c: reconstruction performance vs. computation time')
legend('Block AMM','Hutch h=10','Uniform Sampling')

% visualization: relative error vs runtime (2nd version)
block_sz = 100;
num_group = ceil(k/block_sz);
p_uniform = ones(1,num_group)/num_group;
index_p = cell(1,num_group);
W = cell(1,num_group);
for i = 1:1:num_group
    index_p{i} =  ((i-1)*block_sz+1):(i*block_sz) ; %natural sequential
    W{i} = A(:,index_p{i})*B(index_p{i},:);
end

repeat_time = 10;

s_list_h10 = floor(linspace(100,2000,20));
result_tbl_h10 = zeros(repeat_time,length(s_list_h10));
for rep_num = 1:1:repeat_time
    disp('h10 repeat round = ')
    disp(rep_num)
    counter = 1;
    for s = s_list_h10 
        disp(counter)
        X_h10 = AMM_coarse_hutch(A,B,s,10,index_p);
        result_tbl_h10(rep_num,counter) = norm(X_h10-T,'fro')/norm(T,'fro');
        counter = counter+1;
    end
end

s_list_uni = floor(linspace(100,4000,20));
result_tbl_uni = zeros(repeat_time,length(s_list_uni));
for rep_num = 1:1:repeat_time
    disp('uni repeat round = ')
    disp(rep_num)
    counter = 1;
    for s = s_list_uni 
        disp(counter)
        X_uniform = zeros(size(A,1),size(B,2));
        UniformSample_list = randsample(num_group,s,true,p_uniform);
        for j = 1:1:s
            t_uniform = UniformSample_list(j);           
            X_uniform = X_uniform + W{t_uniform}/p_uniform(t_uniform);                    
        end
        X_uniform = X_uniform/s;
        result_tbl_uni(rep_num,counter) = norm(X_uniform-T,'fro')/norm(T,'fro');
        counter = counter+1;
    end
end

s_list_opt = floor(linspace(100,2000,20));
result_tbl_opt = zeros(repeat_time,length(s_list_opt));
for rep_num = 1:1:repeat_time
    disp('opt repeat round = ')
    disp(rep_num)
    counter = 1;
    for s = s_list_opt 
        disp(counter)
        X_optimal = AMM_true(A,B,s,index_p);
        result_tbl_opt(rep_num,counter) = norm(X_optimal-T,'fro')/norm(T,'fro');
        counter = counter+1;
    end
end

% time measurements
time_tbl1 = zeros(1,length(s_list_opt));
time_tbl2 = zeros(1,length(s_list_h10));
time_tbl3 = zeros(1,length(s_list_uni));

counter = 1;
for s = s_list_opt
    disp(counter)
    f1 = @() AMM_true_tracefun(A,B,s,index_p);
    time_tbl1(counter) = timeit(f1);
    counter = counter+1;
end

counter = 1;
for s = s_list_h10
    disp(counter)
    f2 = @() AMM_coarse_hutch(A,B,s,10,index_p);
    time_tbl2(counter) = timeit(f2);
    counter = counter+1;
end

counter = 1;
for s = s_list_uni
    disp(counter)
    f3 = @() AMM_coarse_uni(A,B,s,index_p);
    time_tbl3(counter) = timeit(f3);
    counter = counter+1;
end

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

save('error_vs_time_bsz100')

% visualization: relative error vs runtime (3rd version)
%block_sz = 100;
block_sz = 500;
num_group = ceil(k/block_sz);
p_uniform = ones(1,num_group)/num_group;
index_p = cell(1,num_group);
%W = cell(1,num_group);
for i = 1:1:num_group
    index_p{i} =  ((i-1)*block_sz+1):(i*block_sz) ; %natural sequential
%    W{i} = A(:,index_p{i})*B(index_p{i},:);
end

%repeat_time = 10;
repeat_time = 5;

% s_list_h10 = floor(linspace(100,2000,20));
s_list_h10 = floor(linspace(10,1000,100));
result_tbl_h10 = zeros(repeat_time,length(s_list_h10));
t_tbl_h10 = zeros(repeat_time,length(s_list_h10));
for rep_num = 1:1:repeat_time
    disp('h10 repeat round = ')
    disp(rep_num)
    counter = 1;
    for s = s_list_h10 
        disp(counter)
        tic
        X_h10 = AMM_coarse_hutch_ver3(A,B,s,10,index_p);
        t_tbl_h10(rep_num,counter)= toc;
        result_tbl_h10(rep_num,counter) = norm(X_h10-T,'fro')/norm(T,'fro');
        counter = counter+1;
    end
end

%s_list_uni = floor(linspace(100,4000,20));
s_list_uni = floor(linspace(10,1000,100));
result_tbl_uni = zeros(repeat_time,length(s_list_uni));
t_tbl_uni = zeros(repeat_time,length(s_list_uni));
for rep_num = 1:1:repeat_time
    disp('uni repeat round = ')
    disp(rep_num)
    counter = 1;
    for s = s_list_uni 
        disp(counter)
        tic
        X_uniform = AMM_coarse_uni_ver2(A,B,s,index_p);
        t_tbl_uni(rep_num,counter)= toc;
%         X_uniform = zeros(size(A,1),size(B,2));
%         UniformSample_list = randsample(num_group,s,true,p_uniform);
%         for j = 1:1:s
%             t_uniform = UniformSample_list(j);           
%             X_uniform = X_uniform + W{t_uniform}/p_uniform(t_uniform);                    
%         end
%         X_uniform = X_uniform/s;
        result_tbl_uni(rep_num,counter) = norm(X_uniform-T,'fro')/norm(T,'fro');        
        counter = counter+1;
    end
end

%s_list_opt = floor(linspace(100,2000,20));
s_list_opt = floor(linspace(10,1000,100));
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

save('error_vs_time_bsz1000_500x1m')