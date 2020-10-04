%% Compare SDA and SE

hold on
[u_inf2, e_inf2, y_inf2, impr2, iteration_number2, error_history2] = SDA(G,d, b,r, u0, R, Q,do_plot);
hold off
xmax =iteration_number2;
xlim([0, xmax])
legend('$e_k$', 'interpreter', 'latex');
set(gca, 'YScale', 'log')
close all 
fig = figure; 
p1 = plotresults(error_history2, iteration_number2, 'SDA', 1, -50, .03, 'Red');
hold on
plot(0:xmax-1, zeros(xmax,1), '--', 'color', 'black');
hold off
legend(p1);