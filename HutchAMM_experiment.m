% Hutch AMM experiment

%% Generating Data
% DATA CASE 1: uniform distribution in the interval (0,1)
m = 50;k = 500;% Define matrix factor size
rng(0)  % For reproducibility
A = rand(m,k);B = rand(k,m);T = A*B;

% DATA CASE 2: Factor A sampled from multivariate normal distribution, with mean of each column linearly spaced increasing.
% Factor B sampled from uniform distribution in the interval (0,1).
mu = 0.1:0.1:50;
Sigma = eye(k);
rng(0) 
A = mvnrnd(mu,Sigma,m);B = rand(k,m);T = A*B;

% DATA CASE 3: Factor A sampled from multivariate normal distribution, with mean of each column sampled from bimodal distribution.
% Factor B sampled from uniform distribution in the interval (0,1).
sigma_a = 1;sigma_b = 2;
mu = [sigma_a*randn([k/2,1]) + 10; sigma_b*randn([k/2,1]) + 30];
Sigma = eye(k);
rng(0)  % For reproducibility
A = mvnrnd(mu,Sigma,m);B = rand(k,m);T = A*B;

%% Initialization & Parameter settings

%%%%%%%%% Define Subgroups/blocks in A and B for AMM

% Subgroups CASE 1: standard AMM - each column forms its subgroup
index_p = cell(k);
for i = 1:1:k
    index_p{i} = i;
end
num_group = length(index_p);

% Subgroups CASE 2: block AMM - subgroup contains multiple multiple
% rows/columns/blocks for A and B
num_group = 50; % set up #subgroups we want
index_p = cell(1,num_group);
for i = 1:1:k
    g_assigned = randsample(num_group,1);
    temp = index_p{g_assigned};
    temp(end+1) = i;
    index_p{g_assigned} = temp; % assign one column index to one group randomly 
end

%%%%%%%%%%%%%% pre-calculate optimal sampling probability for AMM
W = cell(1,num_group); % create a cell array to store A and B subgroups multiplication result
for j = 1:1:num_group
    W{j} = A(:,index_p{j})*B(index_p{j},:);
    p_opt(j) = norm(A(:,index_p{j}))*norm(B(index_p{j},:));
end
p_opt = p_opt/sum(p_opt);

%%%%%%%%%%%%%% Paramters settings
s = 100; % sample number of AMM
num_repeat = 100; % overall experiments repeated times
v = [1,2,5,10,15,20,25,30,40,50]; % #mat-vec multiplication in Hutch estimation
opt_AMM = zeros(1,num_repeat);% for recording relative error for optimal AMM sampling probability
hutch_AMM = zeros(10,num_repeat);% for recording relative error for (length(v) types) Hutch estimated AMM sampling probability


%% Hutch AMM
for i = 1:num_repeat % number of repeated test
    disp(i)
    % calculate AMM estimation with optimal probability
    optIndSample_list = randsample(num_group,s,true,p_opt);
    Tr2 = zeros(size(A,1),size(B,2));
    for j = 1:1:s
        t = optIndSample_list(j);
        Tr2 = Tr2 + W{t}/p_opt(t);
    end
    Tr_Opt = Tr2/s; % optimal AMM estimation
    opt_AMM(i) = norm(Tr_Opt-T,'fro')/norm(T,'fro'); %record its relative error
    % calculate AMM estimation with Hutch estimated probability
    for h = 1:length(v) % loop through diff mat-vec settings for Hutch probability
        % calculate hutch estimator probability
        for j = 1:1:num_group
            p(j) = simple_hutchinson(transpose(W{j})*W{j},v(h));
            p(j) = sqrt(p(j));
        end
        p = p/sum(p);
        % run AMM with hutch probability
        hutchIndSample_list = randsample(num_group,s,true,p);
        Tr = zeros(size(A,1),size(B,2));
        for j = 1:1:s
            t = hutchIndSample_list(j);
            Tr = Tr + W{t}/p(t);
            %Tr = Tr + W{t}/(1/k);
        end
        Tr_Hutch = Tr/s;      
        hutch_AMM(h,i) = norm(Tr_Hutch-T,'fro')/norm(T,'fro');        
    end
end

%% Visualization
figure();
row_names = arrayfun(@num2str,v,'uni',0);
row_names = [row_names 'opt'];
boxplot([hutch_AMM' opt_AMM'],'Labels',row_names)
ylim([0 0.2])
xlabel('number of matrix-vector multiplication')
ylabel('Frobenius relative approximation error')
title('AMM col-group sampling(m = 50,k=500,s=100):Hutch estimation & optimal estimation')
