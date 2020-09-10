%Steepest Decent Algorithm
%Choose K0 = G*
%beta is constant
%Input: G: Plant Model in supervector form
%       d: additional term in Plant Model: y = Gu + d
%       K0: Iterative Matrix
%       beta: design parameter, 0<beta<2/||G||^2
%       u0: start input signal
%       R and Q: positive definite matrices, define inner products in input und
%       output space;length(R) = legnth(u0); length(Q) = length(r); 
%Output: u_inf: final input signal
%        e_inf: final Error, scalar
%        y_inf: final output signal
%        impr: improvement: ||e0||/e_inf, wehere e0 is the initial error
%        iteration_number: number of iterations needed
function [u_inf, e_inf, y_inf, impr,iteration_number] = SDA_beta_as_evs(G,d, r, u0, R, Q,do_plot)


G_star = (R\eye(length(R)))*G'*Q;
G_norm_2 = norm(full(G))^2;

evs = eig(G*G_star); 
evs = evs(evs>10^-12);
evs = evs(evs>.5*G_norm_2);

beta_vec = sort(evs);
beta_vec = 1./beta_vec;
beta_const = 1/G_norm_2;

[u_inf, e_inf, y_inf, impr,iteration_number] = CGA_beta(G,d, G_star,r,beta_vec, beta_const, u0, do_plot);

end

