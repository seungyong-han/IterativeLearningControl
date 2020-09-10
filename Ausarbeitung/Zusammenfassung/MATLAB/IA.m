%% Inverse Model Algorithm
%K0 = G_L, where G_L is the left inverse of the matrix G, if G_L exists, 
%K0 = G_R, where G_R is the left inverse of the matrix G, if G_R exists, 
%Default: left inverse is used 
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
function [u_inf, e_inf, y_inf, impr,iteration_number] = IA(G,d, beta,r, u0, R, Q, do_plot)

if rank(full(G))==length(G(:,1))
    disp('Use left inverse model algorithm')
    [u_inf, e_inf, y_inf, impr,iteration_number] = LIA(G,d, beta,r, u0, R, Q, do_plot)
else
    if rank(full(G))==length(G(1,:))
        [u_inf, e_inf, y_inf, impr,iteration_number]  = RIA(G,d, beta,r, u0, R, Q, do_plot)
        disp('Use right inverse model algorithm')
    else
        disp('Error: Matrix does not have full row or column rank. Can not compute right or left inverse.' )
         u_inf = -1;
         e_inf = -1;
         y_inf = -1;
         impr = -1;
         iteration_number=-1;
         
    end
    
end

end