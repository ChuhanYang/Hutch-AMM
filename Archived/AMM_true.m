function [Tr,m_list] = AMM_true(A,B,m,index_p)
% function AMM (Approximate Matrix Multiplication) uses sampling algorithm.
%  Inputs:
%   m: number of sampling/outer product (usually m << order of T)

k = length(index_p);
Tr = zeros(size(A,1),size(B,2));
% Calculate probability p
p = zeros(k,1);
for j = 1:1:k
    p(j) = norm(A(:,index_p{j})*B(index_p{j},:));
end
p = p/sum(p);
m_list = randsample(k,m,true,p);
for j = 1:1:m
    t = m_list(j);
    Tr = Tr + A(:,index_p{t})*B(index_p{t},:)/p(t);
end
Tr = Tr/m;