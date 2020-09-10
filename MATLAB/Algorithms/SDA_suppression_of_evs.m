%% Steepest descent algorithm with suppression of eigenvalues
%Input: G: Plant Model in supervector form
%       d: additional term in Plant Model: y = Gu + d
%       beta_vec: design parameter for first iterations, vectors are
%       allowed
%       beta_const: design parameter for the next iterations, constant
%       
%        u0: start input signal
%Output: u_inf: cell array, input signal sequence over each 5 iterations
%        e_inf: cell array, includes the error norm over each 5 iterations
%        y_inf: cell array, output signal sequence over each 5 iterations
%        impr: improvement: ||e0||/e_inf{end}, wehere e0 is the initial error
%        iteration_number: number of iterations needed
%        error_history: vector, error norm for each iteration (e_inf is a
%        subset of error_history)
function [u_inf, e_inf, y_inf, impr,iteration_number,error_history] = SDA_suppression_of_evs(G,d, r, u0, R, Q,do_plot)


G_star = sparse((R\eye(length(R)))*G'*Q);
G_norm_2 = norm(full(G*G_star));
beta_vec = linspace(.5*G_norm_2+.5*G_norm_2*10^-10, G_norm_2, 1000);
beta_vec = 1./beta_vec;
beta_const = 1/G_norm_2;
[u_inf, e_inf, y_inf, impr,iteration_number,error_history] = CGA_beta(G,d, G_star,r,beta_vec, beta_const, u0, do_plot);
end

