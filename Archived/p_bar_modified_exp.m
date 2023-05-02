%% Large Matrix Exp with Block Hutch AMM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data generation
% case: uniform
m = 100;k = 10000;% Define matrix factor size
rng(0)  % For reproducibility
A = rand(m,k);B = rand(k,m);T = A*B;

% more meaningful matrices: low rank matrices
m1=100;m2=100;k=1000;k_x=40;k_y=40;%case 1
%m1=1000;m2=1000;k=10000;k_x=400;k_y=400;%case 2
%m1=100;m2=100;k=10000;k_x=40;k_y=40;%case 3
U1 = normrnd(0,1,[m1,k_x]);U2 = normrnd(0,1,[m2,k_y]);
S1 = diag(1-linspace(0,k_x-1,k_x)/k_x);S2 = diag(1-linspace(0,k_y-1,k_y)/k_y);
[temp,~]=qr(randn(k));V1 = temp(1:k_x,:); %unitary matrix
[temp,~]=qr(randn(k));V2 = temp(1:k_y,:);
A = U1*S1*V1;N1=normrnd(0,1,size(A))/m1;
B = (U2*S2*V2)';N2=normrnd(0,1,size(B))/m2;
%A = (A+N1);B = (B+N2);%add noise
%B = A'; %covariance matrix
T = A*B;

% generate exponential decreasing column data
m = 100;k = 10000;
%mu = exp(-(0:0.01:99.99));
mu = exp(linspace(5,0,10000));
Sigma = eye(k);
rng(0) 
A_decay = mvnrnd(mu,Sigma,m);
%A = A_decay;
A = A_decay(:,randperm(k));
B = rand(k,m);T = A*B;

% generate linear increasing
m = 100;k = 10000;
mu = 0:0.001:9.999;
Sigma = eye(k);
rng(0) 
A_decay = mvnrnd(mu,Sigma,m);
%A = A_decay;
A = A_decay(:,randperm(k));
B = 5*rand(k,m);T = A*B;

% generate drift data (exponential decreasing + linear increasing)
m = 100;k = 10000;
mu = 0:0.01:99.99;
mu = exp(-mu);
Sigma = eye(k/2);
rng(0) 
A_decay1 = mvnrnd(mu(1:(k/2)),Sigma,m);
A_decay2 = mvnrnd(linspace(0,1,k/2),Sigma,m);
A_decay = cat(2,A_decay1,A_decay2);
%A = A_decay;
A = A_decay(:,randperm(k));
B = rand(k,m);T = A*B;

% Visualization of column mean
mean_v_1 = zeros(1,k);
mean_v_2 = zeros(1,k);
for i = 1:1:k
    mean_v_1(i) = mean(A_decay(:,i));
    mean_v_2(i) = mean(A(:,i));
end
figure();
subplot(1,2,1);
plot(linspace(1,k,k),mean_v_1)
title('Subplot 1: Before permutation')
subplot(1,2,2);
plot(linspace(1,k,k),mean_v_2)
title('Subplot 2: After permutation')
sgtitle('mean value of each column from A (exponential decreasing case 2)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hutch AMM experiment

% sort the (Hutch estimated) column-wise  (% Enhanced pairwise partition, ditched idea)
W_colw = cell(1,k);
p_Hutch = zeros(1,k);
for j = 1:1:k
    W_colw{j} = A(:,j)*B(j,:);
    p_Hutch(j) = sqrt(simple_hutchinson(transpose(W_colw{j})*W_colw{j},10));
end
[~,idx] = sort(p_Hutch);

% partition based on the sorted or natural sequential
block_sz = 10;
num_group = ceil(k/block_sz);
index_p = cell(1,num_group);
W = cell(1,num_group);
%idx = randperm(k); % random permutate natural index
for i = 1:1:num_group
    %index_p{i} = idx( ((i-1)*block_sz+1):(i*block_sz) ); %sorted results
    index_p{i} =  ((i-1)*block_sz+1):(i*block_sz) ; %natural sequential
    %index_p{i} = idx( ((i-1)*block_sz+1):(i*block_sz) );
    W{i} = A(:,index_p{i})*B(index_p{i},:);
    p_opt(i) = norm(A(:,index_p{i}))*norm(B(index_p{i},:));
end
p_opt_raw = p_opt;
p_opt = p_opt/sum(p_opt);

num_repeat = 100;% number of repeated test
hutch_AMM = zeros(10,num_repeat);
hutch_AMM_bar = zeros(10,num_repeat);
hutch_AMM_plus = zeros(10,num_repeat);
opt_AMM = zeros(1,num_repeat);
v = 3*(1:10);
s = 500; % sample number of AMM
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
    %opt_AMM(i) = norm(Tr_Opt-T,'fro')/norm(T,'fro');
    opt_AMM(i) = norm(Tr_Opt-T,'fro')/(norm(A,'fro')*norm(B,'fro'));
    for h = 1:10 % sample number of Hutch
        % calculate hutch estimator probability
        p = zeros(1,num_group);
        p_plus = zeros(1,num_group);
        for j = 1:1:num_group
            p_plus(j) = hutchplusplus(transpose(W{j})*W{j},v(h));
            p(j) = simple_hutchinson(transpose(W{j})*W{j},v(h));
            p(j) = sqrt(p(j));
            p_plus(j) = sqrt(p_plus(j));
        end
        %p_bar = p_opt_raw/sum(p);
        p_bar_raw = p;
        p = p/sum(p);
        p_plus = p_plus/sum(p_plus);
        % run AMM with hutch probability
        hutchIndSample_list = randsample(num_group,s,true,p);
        hutchPlusIndSample_list = randsample(num_group,s,true,p_plus);
        Tr = zeros(size(A,1),size(B,2));
        Tr_bar = zeros(size(A,1),size(B,2));
        Tr_plus = zeros(size(A,1),size(B,2));
        sampledset = zeros(1,s);
        for j = 1:1:s
            t = hutchIndSample_list(j);
            t_plus = hutchPlusIndSample_list(j);            
            t_bar = randsample(num_group,1,true,p_bar_raw/sum(p_bar_raw));                        
           
            Tr = Tr + W{t}/p(t);
            Tr_plus = Tr_plus + W{t_plus}/p_plus(t_plus);
            % replace the estimated norm with sampled true norm
            if ismember(t_bar,sampledset) ~= true
                p_bar_raw(t_bar) = norm(W{t_bar},'fro');
            end
            sampledset(j) = t_bar;
            Tr_bar = Tr_bar + W{t_bar}/(p_bar_raw(t_bar)/sum(p_bar_raw));                     
        end
        Tr_Hutch = Tr/s;
        Tr_Hutch_bar = Tr_bar/s;
        Tr_Hutch_plus = Tr_plus/s; 
        %hutch_AMM(h,i) = norm(Tr_Hutch-T,'fro')/norm(T,'fro');
        hutch_AMM(h,i) = norm(Tr_Hutch-T,'fro')/(norm(A,'fro')*norm(B,'fro'));
        %hutch_AMM_bar(h,i) = norm(Tr_Hutch_bar-T,'fro')/norm(T,'fro');
        hutch_AMM_bar(h,i) = norm(Tr_Hutch_bar-T,'fro')/(norm(A,'fro')*norm(B,'fro'));
        %hutch_AMM_plus(h,i) = norm(Tr_Hutch_plus-T,'fro')/norm(T,'fro');
        hutch_AMM_plus(h,i) = norm(Tr_Hutch_plus-T,'fro')/(norm(A,'fro')*norm(B,'fro'));
    end
end

% visualization of three methods comparison
figure();
row_names = arrayfun(@num2str,v,'uni',0);
row_names = [row_names 'opt'];
subplot(3,1,1);
boxplot([hutch_AMM' opt_AMM'],'Labels',row_names)
%ylim([0.011 0.018])
title('Subplot 1: p_{Hutch}')
subplot(3,1,2);
boxplot([hutch_AMM_bar' opt_AMM'],'Labels',row_names)
%ylim([0.011 0.018])
title('Subplot 2: p_{bar}')
subplot(3,1,3);
boxplot([hutch_AMM_plus' opt_AMM'],'Labels',row_names)
%ylim([0.011 0.018])
xlabel('number of matrix-vector multiplication')
title('Subplot 3: p_{Hutch++}')
sgtitle('exponential decreasing case 2, sample number = 500, block size = 10')

% check if the code is correct
m = 100;k = 1000;rng(0);A = rand(m,k);B = rand(k,m);T = A*B;
block_sz = 1; % fix block size = 1
num_group = ceil(k/block_sz);
index_p = cell(1,num_group);
W = cell(1,num_group);
%idx = randperm(k); % random permutate natural index
for i = 1:1:num_group
    index_p{i} =  ((i-1)*block_sz+1):(i*block_sz) ; %natural sequential
    %index_p{i} = idx( ((i-1)*block_sz+1):(i*block_sz) );
    W{i} = A(:,index_p{i})*B(index_p{i},:);
    p_opt(i) = norm(A(:,index_p{i}))*norm(B(index_p{i},:));
end
p_opt_raw = p_opt;
p_opt = p_opt/sum(p_opt);

appr_tbl = zeros(20,100);
opt_AMM = zeros(1,100);
counter = 1;
for s = linspace(10,1000,100)% sample number of AMM
    disp(counter)
    for h = 1:20 % number of repeat
        % calculate hutch estimator probability
        p = zeros(1,num_group);
        for j = 1:1:num_group
            p(j) = simple_hutchinson(transpose(W{j})*W{j},1);
            p(j) = sqrt(p(j));
        end
        p = p/sum(p);
        % run AMM with hutch probability
        hutchIndSample_list = randsample(num_group,s,true,p);
        Tr = zeros(size(A,1),size(B,2));
        for j = 1:1:s
            t = hutchIndSample_list(j);                          
            Tr = Tr + W{t}/p(t);
        end
        Tr_Hutch = Tr/s;
        %appr_tbl(h,counter) = norm(Tr_Hutch-T,'fro')/norm(T,'fro');
        appr_tbl(h,counter) = norm(Tr_Hutch-T,'fro')/(norm(A,'fro')*norm(B,'fro'));
    end
    % optimal case
    optIndSample_list = randsample(num_group,s,true,p_opt);
    Tr_opt = zeros(size(A,1),size(B,2));
    for j = 1:1:s
       t = optIndSample_list(j);                          
       Tr_opt = Tr_opt + W{t}/p_opt(t);
    end
    Tr_opt = Tr_opt/s;
    %opt_AMM(1,counter) = norm(Tr_opt-T,'fro')/norm(T,'fro');
    opt_AMM(1,counter) = norm(Tr_opt-T,'fro')/(norm(A,'fro')*norm(B,'fro'));
    
    counter = counter + 1;
end

figure();
row_names = arrayfun(@num2str,linspace(10,1000,100),'uni',0);
boxplot(appr_tbl,'Labels',row_names)
%ylim([0 0.5])
ylim([0 0.8])
hold on
plot(1:1:100,opt_AMM)
title('Increasing AMM sampling number with fixed box size = 1, low rank case')
hold off