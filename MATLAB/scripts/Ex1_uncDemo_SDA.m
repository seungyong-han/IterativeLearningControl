%% Uncertainty Demonstration 

A = [2 1; 4 3];
[G, d] = get_G(A - B*F, B, C - D*F, D, x0, N-1);
K0 = sparse((R\eye(length(R)))*G'*Q);
beta = b*2/max(abs(eig(G*K0)));


A = [2 1; 3 3];
[G, d] = get_G(A - B*F, B, C - D*F, D, x0, N-1);

[u_inf1, e_inf1, y_inf1, impr1, iteration_number1, error_history1] = CGA(G,d, K0,beta,r, u0, do_plot);
close all
hold on 
plot(0:N-1, r);
y = y_inf1{end};
plot(0:N-1, y);
legend('reference signal', 'solution');
hold off
