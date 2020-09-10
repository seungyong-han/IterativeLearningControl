%% Right Inverse Model Algorithm
%K0 = G_L, where G_L is the left inverse of the matrix G 
%Here: Moore-Penrose pseudoinverse is used
%beta=const
%Input: G: Plant Model in supervector form
%       d: additional term in Plant Model: y = Gu + d
%       K0: Iterative Matrix
%       beta: design parameter, 0<beta<1; Here: beta is the quota of max.
%       possible beta to be chosen. For beta = 1/2 get solution in one
%       iteration, but very bad robustness properties. 
%       u0: start input signal
%       R and Q: positive definite matrices, define inner products in input und
%       output space;length(R) = legnth(u0); length(Q) = length(r); 
%Output: u_inf: final input signal
%        e_inf: norm of the final error
%        y_inf: final output signal
%        impr: improvement: ||e0||/e_inf, wehere e0 is the initial error
%        iteration_number: number of iterations needed
function [u_inf, e_inf, y_inf, impr,iteration_number] = LIA(G,d, beta,r, u0, R, Q, do_plot)

if rank(full(G))~=length(G(:,1))
    disp('Error: G has not full column rank. Can not compute left inverse')
    u_inf = -1;
    e_inf = -1;
    y_inf = -1;
    impr = -1;
    iteration_number = -1;
else
    K0 = pinv(full(G)')';
    K0 = sparse(K0);
    beta_max = 2/norm(full(G))^2;
    beta = beta*beta_max
    [u_inf, e_inf, y_inf, impr,iteration_number] = CGA(G,d, K0,beta,r, u0, do_plot);
end

end

