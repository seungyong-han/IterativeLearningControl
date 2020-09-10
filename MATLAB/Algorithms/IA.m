%% Inverse Model Algorithm
%K0 = G_R, where G_R is the left inverse of the matrix G, if G_R exists, 
%K0 = G_L, where G_L is the left inverse of the matrix G, if G_L exists, 
%Default: right inverse is used 
%Here: use Moore-Penrose pseudoinverse
%beta=const
%Input: G: Plant Model in supervector form
%       d: additional term in Plant Model: y = Gu + d
%       K0: Learning Matrix
%       beta: design parameter, 0<beta<1; Here: beta is the quota of max.
%       possible beta to be chosen. For beta = 1/2 get solution in one
%       iteration, but bad robustness properties. 
%       u0: initial input signal
%Output: u_inf: cell array, input signal sequence over each 5 iterations
%        e_inf: cell array, includes the error norm over each 5 iterations
%        y_inf: cell array, output signal sequence over each 5 iterations
%        impr: improvement: ||e0||/e_inf{end}, wehere e0 is the initial error
%        iteration_number: number of iterations needed
%        error_history: vector, error norm for each iteration (e_inf is a
%        subset of error_history)
function [u_inf, e_inf, y_inf, impr,iteration_number, error_history] = IA(G,d, beta,r, u0, do_plot, do_print)



if rank(full(G))==length(G(:,1))
    if(do_print)
        [u_inf, e_inf, y_inf, impr,iteration_number, error_history]  = RIA(G,d, beta,r, u0, do_plot);
        disp('Use right inverse model algorithm')
    else
        [u_inf, e_inf, y_inf, impr,iteration_number, error_history]  = RIA(G,d, beta,r, u0, do_plot);
    end
elseif rank(full(G))==length(G(1,:))
    if(do_print)
        disp('Use left inverse model algorithm')
        [u_inf, e_inf, y_inf, impr,iteration_number, error_history] = LIA(G,d, beta,r, u0,  do_plot);
    else
        [u_inf, e_inf, y_inf, impr,iteration_number, error_history] = LIA(G,d, beta,r, u0, do_plot);
    end
else
        fprintf(2, 'Error: Matrix G does not have full row or column rank. Can not compute right or left inverse. \n' )
         u_inf = -1;
         e_inf = -1;
         y_inf = -1;
         impr = -1;
         iteration_number=-1;
         error_history = -1; 
end
    
end
