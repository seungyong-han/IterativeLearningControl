%% IA vs LQR

close all
[u_inf1, e_inf1, y_inf1, impr1, iteration_number1, error_history1] =  RIA(G,d, b,r_vec, u0, 0);%IA(G,d, b,r, u0, do_plot, 1);
xmax = iteration_number1;
if(xmax>0)
    xlim([0, xmax])
end
legend('$e_k$', 'interpreter', 'latex');

%% Plot results LQR vs IA 
close all
plot(0:N, r_vec);
hold on
plot(0:N, y_sv);
plot(0:N, y_inf1{end}, '--');
legend('Reference','LQR', 'PIA');
%legend('Reference', 'IA');
xlim([0, N-1])
hold off
