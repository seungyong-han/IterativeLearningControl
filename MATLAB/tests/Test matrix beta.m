%% Test matrix beta

clc
clear all


rand_sys = drss(100, 5, 7);
rand_sys = prescale(rand_sys);
rand_sys = minreal(rand_sys);
[A, B, C, D, ts] = ssdata(rand_sys);


%rank(D, Dk) must be min{m,l}
x0 = zeros(length(A), 1);
N = 10;
order = 40;


%m-output l-input system; in supervector description m*(N+1) - output
%l*(N+1)-input
l  = length(B(1,:));
m = length(C(:,1));


u0 = repmat(0.1,l*(N+1),1);
%r = G*G'*5*rand(12,1) + G*u0 + d;%to ensure that e0 is in image(GG*)
r = 5*rand(m*(N+1),1);

R = eye(l*(N+1));

Q = eye(m*(N+1));
[G, d] = get_G(A, B, C, D, x0, N);


beta = .9;

do_plot = 1;




[G, d] = get_G(A, B, C, D, x0, N);
[Ak, Bk, Ck, Dk, x0k] = get_reduced_system(A, B, C, D, x0, order, 1);

K = get_G(Ak, Bk, Ck, Dk, x0k,  N);
K0 = get_adj(R, K, Q);
G_star = get_adj(R, G, Q);

beta = beta*2/norm(full(K0*G));




%if rank(full(G*K'))<req_rankK
if sum(eig(full(G*K0 + K*G_star - beta*K*G_star*G*K0))<0)>0 %Monotonic convergence is guaranteed if GK* + KG* > beta KG*GK* 
    disp('Error: ker(GK*) neq {0}')
    u_inf = -1;
    e_inf = -1;
    y_inf = -1;
    impr = -1;
    iteration_number = -1;
else
    
    K_old = get_adj(R, K, Q);
    K0 = K_old; 
    M = full(G*K0);
    beta_1 = ones(10, 1)*1/norm(M(:,1:10));
    beta_2 = ones(dim_M-20,1)*1/norm(M(:,11:dim_M-10));
    beta_3 = ones(10,1)*1/norm(M(:,dim_M-9:dim_M));
    
    beta = diag([beta_1; beta_2; beta_3]);
        
%     beta = ones(length(M),1);
%     for i = 1:length(M)
%         beta(i) = 1/norm(M(i,:));
%     end
%     beta = diag(beta);
    
     
     K0 = K0*beta;
     M = G*K0;
     beta = 1/norm(M);
    
 
%%    
 
     beta = 1;
     e = e0;
     u = u0;
     betaK0 = beta*K0;
     betaGK0 = beta*G*K0;
     while cont
        iteration_number = iteration_number + 1;
        u_new = u + betaK0*e;
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