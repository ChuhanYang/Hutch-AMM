% experiments on natural data
unit_change = 1000/3600;
A = X(:,1:2700)*unit_change;
rou_jam = 1/7.5;
c_cong = -15*unit_change;
c_free = 100*unit_change;
B = real(rou_jam * (1 + c_free/c_cong * log(1-A/c_free)).^(-1))';
[m,k] = size(A);T=A*B;

% TF-IDF TDT2 data
load("TDT2.mat")
fea_sample = fea(randsample(size(fea,1),1000,false),:);
fea = tfidf(fea_sample);
A = full(fea);[m,k] = size(A);B = A';
T=A*B;

% Movie property data
load("X_movie.mat")
A = double(count_matrix);[m,k] = size(A);B=A';
T=A*B;

% Wiki People data
load("X_wikipeople.mat")
A = full(wiki_matrix);A = A(randsample(size(A,1),1000,false),:);[m,k] = size(A);B=A';
T=A*B;

% 20news data
load("X_20news.mat")
% readtable('X_20news.mat')
% dlmread("X_20news.mat")
A = full(newsgroup);A = A(randsample(size(A,1),1000,false),:);[m,k] = size(A);B=A';
T=A*B;

% recv1 data
load("recv1.mat")
A = recv1_matrix;A = full(A(randsample(size(A,1),1000,false),:));[m,k] = size(A);B=A';
T=A*B;

% Beijing weather data
A = readtable('PRSA_Data_Wanshouxigong_20130301-20170228.csv');
A = PRSADataChangping2013030120170228(:,[6:15,17]);
A = table2array(A)';
A = fillmissing(A,'previous',2);
%A = normr(A);
[m,k] = size(A);B=A';T=A*B;

% household electric power consumption data
A = householdpowerconsumption(:,3:end);
A = table2array(A)';
A = fillmissing(A,'previous',2);
[m,k] = size(A);B=A';T=A*B;

% yellow taxi data
A = yellowtaxinum(:,2:end);
A = table2array(A)';
[m,k] = size(A);B=A';T=A*B;

%block_sz = 100;
%block_sz = 500;
%block_sz = 24;
%block_sz = 60;
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

% s_list_h10 = floor(linspace(10,300,10)); % TF-IDF
% s_list_h10 = floor(linspace(10,1000,10)); % weather
% s_list_h10 = floor(linspace(10,1000,5)); % synthetic MA
% s_list_h10 = floor(linspace(10,7000,5));
% s_list_h10 = floor(linspace(10,20000,5)); % yellow taxi (old)
s_list_h10 = floor(linspace(10,30000,5));
%s_list_h10 = floor(linspace(1000,400000,10));
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


% s_list_uni = floor(linspace(10,300,10)); % TF-IDF
% s_list_uni = floor(linspace(10,1000,10)); % weather
% s_list_uni = floor(linspace(10,1000,5)); % synthetic MA
% s_list_uni = floor(linspace(10,7000,5));
s_list_uni = floor(linspace(10,60000,5)); % yellow taxi
%s_list_uni = floor(linspace(1000,400000,10));
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
s_list_opt = floor(linspace(10,30000,5)); % yellow taxi
%s_list_opt = floor(linspace(1000,400000,10));
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

% % non-coarse optimal AMM
% s_list_AMM = floor(linspace(10,0.2*500*block_sz,10));
% result_tbl_AMM = zeros(repeat_time,length(s_list_AMM));
% t_tbl_AMM = zeros(repeat_time,length(s_list_AMM));
% for rep_num = 1:1:repeat_time
%     disp('standard AMM repeat round = ')
%     disp(rep_num)
%     counter = 1;
%     for s = s_list_AMM 
%         disp(counter)
%         tic
%         X_optimal = AMM_standard_fair(A,B,s,block_sz);
%         t_tbl_AMM(rep_num,counter)= toc;
%         result_tbl_AMM(rep_num,counter) = norm(X_optimal-T,'fro')/norm(T,'fro');
%         counter = counter+1;
%     end
% end

% non-coarse optimal AMM
% block_sz = 1;
% %block_sz = 500;
% num_group = ceil(k/block_sz);
% p_uniform = ones(1,num_group)/num_group;
% index_p = cell(1,num_group);
% %W = cell(1,num_group);
% for i = 1:1:num_group
%     index_p{i} =  ((i-1)*block_sz+1):(i*block_sz) ; %natural sequential
% %    W{i} = A(:,index_p{i})*B(index_p{i},:);
% end
% 
% %s_list_AMM = floor(linspace(10,2*300*block_sz,10));  % TF-IDF
% %s_list_AMM = floor(linspace(10,2*150*block_sz,10)); % newsgroup
% s_list_AMM = floor(linspace(10,2*500*block_sz,10));
% result_tbl_AMM = zeros(repeat_time,length(s_list_AMM));
% t_tbl_AMM = zeros(repeat_time,length(s_list_AMM));
% for rep_num = 1:1:repeat_time
%     disp('standard AMM repeat round = ')
%     disp(rep_num)
%     counter = 1;
%     for s = s_list_AMM 
%         disp(counter)
%         tic
%         X_optimal = AMM_true_tracefun_ver2(A,B,s,index_p);
%         t_tbl_AMM(rep_num,counter)= toc;
%         result_tbl_AMM(rep_num,counter) = norm(X_optimal-T,'fro')/norm(T,'fro');
%         counter = counter+1;
%     end
% end

% % result - time of standard matrix multiple
% std_time_tbl = zeros(1,repeat_time);
% for rep_num = 1:1:repeat_time
%     tic
%     T = A*B;
%     std_time_tbl(rep_num) = toc;
%     clear T;
% end
% std_time = mean(std_time_tbl);


