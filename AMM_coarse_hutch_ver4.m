function [Tr,m_list] = AMM_coarse_hutch_ver4(A,B,m,c,index_p)
% function AMM (Approximate Matrix Multiplication) but uses hucthplusplus/simple_hutchinson
% to estimate sampling probability. AMM_toep_coarse Samples column groups/corase partion 
% instead of fineset partition 
% Considers subgroup with diff size
%  Inputs:
%   A,B: target matrix AB should be estimated
%   m: number of sampling/outer product in AMM
%   c: hyperparameter from hucthplusplus/simple_hutchinson
%   index_p: cell arrays of index partition, eg. {[1,3,5],2,[2]}

k = length(index_p);
Tr = zeros(size(A,1),size(B,2));
% Calculate probability p
p = zeros(k,1);
% for j = 1:1:k
%     W = A(:,index_p{j})*B(index_p{j},:);
%     p(j) = sqrt(simple_hutchinson(transpose(W)*W,c));
% end
if length(index_p{1})<size(A,1) % change matrix multiplication chain based on whether the matrix is wide or thin.
    S = 2*randi(2,length(index_p{1}),c)-3;
    for j = 1:1:(k-1)
        p(j) = sqrt(trace(S'*A(:,index_p{j})'*A(:,index_p{j})*B(index_p{j},:)*B(index_p{j},:)'*S) / c);
    end
    S_last = 2*randi(2,length(index_p{k}),c)-3;
    p(k) = sqrt(trace(S_last'*A(:,index_p{k})'*A(:,index_p{k})*B(index_p{k},:)*B(index_p{k},:)'*S_last) / c);
else
   S = 2*randi(2,size(A,1),c)-3;
   for j = 1:1:k       
        p(j) = sqrt(trace(S'*A(:,index_p{j})*B(index_p{j},:)*B(index_p{j},:)'*A(:,index_p{j})'*S) / c);
   end
end

p = real(p/sum(p));
m_list = randsample(k,m,true,p);
for j = 1:1:m
    t = m_list(j);
    Tr = Tr + (A(:,index_p{t})*B(index_p{t},:))/p(t);
end
Tr = Tr/m;