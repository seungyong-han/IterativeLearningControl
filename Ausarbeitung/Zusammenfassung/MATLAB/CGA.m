%%
%Conceptual Gradient Algorithm for constant beta -> Can also be used for
%inverse algorithms 
%Input: G: Plant Model in supervector form
%       d: additional term in Plant Model: y = Gu + d
%       K0: Iterative Matrix
%       beta: design parameter, 0<beta<2/||G||^2
%       u0: start input signal
%Output: u_inf: final input signal
%        e_inf: final Error, scalar
%        y_inf: final output signal
%        impr: improvement: ||e0||/e_inf, wehere e0 is the initial error
%        iteration_number: number of iterations needed
function [u_inf, e_inf, y_inf, impr, iteration_number] = CGA(G,d, K0,beta,r, u0, do_plot)
u = u0;
e0 = r - G*u - d; 
e = e0;
cont = 1;
iteration_number = 0;
error_history = [norm(e0)];
input_history = [norm(u0)];

betaGK0 = beta*G*K0;%calculate once
while cont
    iteration_number = iteration_number + 1;
    u_new = u + beta*K0*e;
    e_new = (eye(length(e0)) - betaGK0)*e;
    
    if norm(e - e_new)<10^-10
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
    
end

u_inf = u_new;
e_inf = norm(e_new);
y_inf = G*u_new + d;
impr = norm(e0)/e_inf;

if do_plot
    plot(0:iteration_number, error_history);
end

end

