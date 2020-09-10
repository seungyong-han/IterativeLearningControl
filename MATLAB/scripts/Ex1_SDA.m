%% Compare SDA and SE

hold on
[u_inf2, e_inf2, y_inf2, impr2, iteration_number2, error_history2] = SDA(G,d, b,r, u0, R, Q,do_plot);
[u_inf3, e_inf3, y_inf3, impr3, iteration_number3, error_history3] = SDA_suppression_of_evs(G,d, r, u0, R, Q,do_plot);
hold off
xmax = max(iteration_number2, iteration_number3);
xlim([0, xmax])
legend('$e_k$', 'interpreter', 'latex');
set(gca, 'YScale', 'log')
close all 
plotresults(error_history2, iteration_number2, 'SDA', 1, -70, 1e-3)
plotresults(error_history3, iteration_number3, 'SDA_ES', 1, 10, 0)

