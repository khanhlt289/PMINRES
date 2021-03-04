% function [x, converged, iter_cnt, res_norms] = PMINRES_repair(A, b, res_tol, max_iter,M)
% Thu nghiem PMINRES algorithm
% Ngay sua chua: 22/2/2021
% Nguoi sua chua: Luu Truong Khanh
% Correspond to Algorithm 53 in Holger Wendland Numerical Linear Algebra An Introduction"
% Pass exam with M = eye(n)! 
% Ma tran tien dieu kien PMINRES co dang doi xung co the la Cholesky

% n = size(A,1);
% l = zeros(n,n);
% l(1,1) = sqrt(A(1,1)); j = 1;
% for i = j + 1 : n
%     l(i,j) = A(i,j);
%     for k = 1:j-1
%         
%         l(i,j) = l(i,j) - l(i,k) * l(j,k);
%     end
%     l(i,j) = l(i,j) / l(j,j);
% end
% for j = 2:n
%     l(j,j) = A(j,j);
%     for k = 1 : j-1
%         l(j,j) = l(j,j) - l(j,k)^2;
%     end
%     l(j,j) = sqrt(l(j,j));
%     for i = j + 1 : n
%         l(i,j) = A(i,j);
%         for k = 1:j-1
%             l(i,j) = l(i,j) - l(i,k) * l(j,k);
%         end
%         l(i,j) = l(i,j) / l(j,j);
%     end
% end

    clear all; clc;
    A = xlsread('D:/Study/FEM Matlab KLT/PTHH/kk_UBoot.xlsx');
    b = xlsread('D:/Study/FEM Matlab KLT/PTHH/ff_UBoot.xlsx');

    n = size(A,1);
    [a_val, a_row_ptr, a_col_idx] = sparse2csr(A,1);
%     b = rand(n,1);
    res_tol = 1e-10;
    max_iter = n;
    
    l = Chol_trial(A);
    % Ma tran tien dieu kien Cholesky cua 
    % he phuong trinh dai so tuyen tinh Ax = b;
    M = l * l';
    %     M = eye(n,n);
    m_val = ones(n,1);
    m_row_ptr = [1:n+1]';
    m_col_idx = [1:n]';
    
%     [m_val, m_row_ptr, m_col_idx] = sparse2csr(M,1);
    m_diagonal = ones(n,1);
    
%     M = tril(A);
%     M = diag(diag(A));
%     M = rand(64,64);
    
%     n = size(A, 1);

%     if (nargin < 4) 
%         res_tol  = 1e-10; end
%     
%     if (nargin < 5) 
%         max_iter = 1000; end
%     
% 	v = cell(1,3);
%     for i = 1:3
%         v{1,i} = zeros(n,1);
%     end
    
	p = cell(1,3);
    for i = 1:3
        p{1,i} = zeros(n,1);
    end
    
    z = cell(1,2);
    for i = 1:2
        z{1,i} = zeros(n,1);
    end
    
    x = zeros(n,1);
    w = zeros(n,1);
    
    beta = zeros(2,1);
    s = zeros(3,1);
    c = ones(3,1);
    
%     r = b - A * x; 
%     r = b - dotproduct(a_val, a_row_ptr, a_col_idx, x); 
%     delta = norm(r); 
%     v{1,2} = r / delta;
    
    r = b - dotproduct(a_val, a_row_ptr, a_col_idx, x); 
    delta = norm(r); 
    v{1,2} = r / delta;
%     z{1,1} = M * v{1,2};
    z{1,1} = dotproduct(m_val, m_row_ptr, m_col_idx, v{1,2});
%     z{1,1} = M * v{1,2};
    gamma = zeros(4,1);
    
    stop_res = delta * res_tol;
    iter_cnt = 1;
    res_norms = zeros(max_iter, 1);
    converged = 0;
%         w = A * v{1,2};
        w = dotproduct(a_val, a_row_ptr, a_col_idx, v{1,2});
        alpha = dot(w, z{1,1});
        w = w - alpha * v{1,2};
%         w_wilde = M * w;
        w_wilde = dotproduct(m_val, m_row_ptr, m_col_idx, w);
        beta(2,1) = sqrt(dot(w, w_wilde));  

        if (abs(delta) > stop_res)
            gamma(1,1) = c(2,1) * alpha - s(2,1) * c(1,1) * beta(1,1);
            gamma(2,1) = (gamma(1,1)^2 + beta(2,1)^2)^0.5;
            gamma(3,1) = s(2,1) * alpha + c(2,1) * c(1,1) * beta(1,1);
            gamma(4,1) = s(1,1) * beta(1,1);
            c(3,1) = gamma(1,1)/ gamma(2,1); 
            s(3,1) = beta(2,1) / gamma(2,1);
            p{1,3} = z{1,1} - gamma(3,1) * p{1,2} - gamma(4,1) * p{1,1};
            p{1,3} = p{1,3} / gamma(2,1);
            x = x + c(3,1) * delta * p{1,3};
            delta = - s(3,1) * delta;
            v{1,3} = w / beta(2,1);
            z{1,2} = w_wilde / beta(2,1);
        end
        res_norms(iter_cnt) = abs(delta);
        iter_cnt = iter_cnt + 1;
    while (iter_cnt < max_iter)
        beta(1,1) = beta(2,1);
        s(1,1) = s(2,1); s(2,1) = s(3,1);
        c(1,1) = c(2,1); c(2,1) = c(3,1);
        p{1,1} = p{1,2}; p{1,2} = p{1,3};
        v{1,1} = v{1,2}; v{1,2} = v{1,3}; 
        z{1,1} = z{1,2};
        w = dotproduct(a_val, a_row_ptr, a_col_idx, v{1,2}) - beta(1,1) * v{1,1};
        alpha = dot(w, z{1,1});
        w = w - alpha * v{1,2};
        w_wilde = dotproduct(m_val, m_row_ptr, m_col_idx, w);
        beta(2,1) = sqrt(dot(w, w_wilde));  

        if (abs(delta) > stop_res)
            gamma(1,1) = c(2,1) * alpha - s(2,1) * c(1,1) * beta(1,1);
            gamma(2,1) = (gamma(1,1)^2 + beta(2,1)^2)^0.5;
            gamma(3,1) = s(2,1) * alpha + c(2,1) * c(1,1) * beta(1,1);
            gamma(4,1) = s(1,1) * beta(1,1);
            c(3,1) = gamma(1,1)/ gamma(2,1); 
            s(3,1) = beta(2,1) / gamma(2,1);
            p{1,3} = z{1,1} - gamma(3,1) * p{1,2} - gamma(4,1) * p{1,1};
            p{1,3} = p{1,3} / gamma(2,1);
            x = x + c(3,1) * delta * p{1,3};
            delta = - s(3,1) * delta;
            v{1,3} = w / beta(2,1);
            z{1,2} = w_wilde / beta(2,1);
        else
            break;
        end
        res_norms(iter_cnt) = abs(delta);
        iter_cnt = iter_cnt + 1;
    end
    if (abs(delta) <= res_tol) 
        converged = 1; 
    end
    	
    res_norms = res_norms(1:iter_cnt);