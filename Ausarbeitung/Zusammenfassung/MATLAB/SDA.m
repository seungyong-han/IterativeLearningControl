%Steepest Decent Algorithm
%Choose K0 = G*
%beta is constant
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
%        e_inf: final Error, scalar
%        y_inf: final output signal
%        impr: improvement: ||e0||/e_inf, wehere e0 is the initial error
%        iteration_number: number of iterations needed
function [u_inf, e_inf, y_inf, impr,iteration_number] = SDA(G,d, beta,r, u0, R, Q, do_plot)


G_star = sparse((R\eye(length(R)))*G'*Q)

beta = beta*2/norm(full(G))^2;
[u_inf, e_inf, y_inf, impr,iteration_number] = CGA(G,d, G_star,beta,r, u0, do_plot);
end

