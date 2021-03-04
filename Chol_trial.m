% Thu nghiem Cholesky bang lenh chol
% tic;
% A = delsq(numgrid('L',1000));  n = 100;
% L = ichol(A,struct('type','ict','droptol',1e-02,'michol','on'));
% e = ones(size(A,2),1);t
% norm(A*e - L*L'*e, 'fro')
% toc;

clear all; clc;
% function l = Chol_trial(A)
A = delsq(numgrid('S',28));
n = size(A,1);
l = zeros(n,n);
l(1,1) = sqrt(A(1,1)); j = 1;
for i = j + 1 : n
    l(i,j) = A(i,j);
    for k = 1:j-1
        l(i,j) = l(i,j) - l(i,k) * l(j,k);
    end
    l(i,j) = l(i,j) / l(j,j);
end

for j = 2:n
    l(j,j) = A(j,j);
    for k = 1 : j-1
        l(j,j) = l(j,j) - l(j,k)^2;
    end
    l(j,j) = sqrt(l(j,j));
    for i = j + 1 : n
        l(i,j) = A(i,j);
        for k = 1:j-1
            l(i,j) = l(i,j) - l(i,k) * l(j,k);
        end
        l(i,j) = l(i,j) / l(j,j);
    end
end