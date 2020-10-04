%% Right Inverse Model Algorithm
%K0 = G_R, where G_R is the right inverse of the matrix G 
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
function [u_inf, e_inf, y_inf, impr,iteration_number, error_history] = RIA(G,d, beta,r, u0, do_plot)


%check if the matrix G has full row rank
% if rank(full(G))~=length(G(:,1))
%     disp('Error: G has not full row rank. Can not compute right inverse')
%     u_inf = -1;
%     e_inf = -1;
%     y_inf = -1;
%     impr = -1;
%     iteration_number = -1;
%     error_history = -1;
%else %write the algorithm down because of other error calculation 
    beta_max = 2;
    beta = beta*beta_max;
    beta = beta*2;
    %K0 = pinv(full(G));
    K0 = G\eye(length(G));
    K0 = sparse(K0);

    u = u0;
    e0 = r - G*u - d; 
    e = e0;
    cont = 1;
    iteration_number = 0;
    error_history = [norm(e0)];
    input_history = [norm(u0)];

    u_inf = {};
    e_inf = {};
    y_inf = {};
    cell_nb = 0;
    while cont
        iteration_number = iteration_number + 1;
        u_new = u + beta*K0*e;
        %e_new = (1 - beta)*e;
        e_new = r - G*u_new - d; %for demonstration non-stability

        if norm(e - e_new)<10^-6
            cont = 0;
        end
        error_history = [error_history, norm(e_new)];
        input_history = [input_history, norm(u_new)];
        if(mod(iteration_number, 10000) == 0)
            disp("currErrordiff:")
            norm(e_new - e)
        end
        u = u_new;
        e = e_new;

        if mod(iteration_number,5) == 0
            cell_nb = cell_nb + 1;
            u_inf{cell_nb} = u_new;
            e_inf{cell_nb} = norm(e_new);
            y_inf{cell_nb} = G*u_new + d;
        end


 %   end

    if mod(iteration_number,5) ~= 0
            cell_nb = cell_nb + 1;
            u_inf{cell_nb} = u_new;
            e_inf{cell_nb} = norm(e_new);
            y_inf{cell_nb} = G*u_new + d;
    end

    impr = norm(e0)/e_inf{length(e_inf)};

    if do_plot
        plot(0:iteration_number, error_history);
    end

end
end



