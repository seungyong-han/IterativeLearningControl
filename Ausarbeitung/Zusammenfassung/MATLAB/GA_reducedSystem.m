%necessary cindition for monotonic convergence: GK* + KG* > beta KG*GK*, in particular ker(G*K0) = {0} 
%INPUT: beta: value between 0 and 1, Percentage of maximum value for coeficient beta

function [u_inf, e_inf, y_inf, impr, iteration_number] = GA_reducedSystem(A, B, C, D, x0, N, order, reducer, beta, r, u0, R, Q,do_plot)

[G, d] = get_G(A, B, C, D, x0, N);
[Ak, Bk, Ck, Dk, x0k] = get_reduced_system(A, B, C, D, x0, order, reducer);

K = get_G(Ak, Bk, Ck, Dk, x0k,  N);
K0 = get_adj(R, K, Q);
G_star = get_adj(R, G, Q);

beta = beta*2/norm(full(K0*G));
%beta = beta/norm(full(G))^2; %beta from p.187




%if rank(full(G*K'))<req_rank(K)
if sum(eig(full(G*K0 + K*G_star - beta*K*G_star*G*K0))<0)>0 %Monotonic convergence is guaranteed if GK* + KG* > beta KG*GK* 
    disp('Error: ker(GK*) neq {0}')
    u_inf = -1;
    e_inf = -1;
    y_inf = -1;
    impr = -1;
    iteration_number = -1;
else
eig((eye(length(r)) - beta*G*K0))%if the evs of GK0 lie in large range => max(eig(L))is almost 1 => slow convergence 
[u_inf, e_inf, y_inf, impr, iteration_number] = CGA(G,d, K0,beta,r, u0, do_plot);
end
end

