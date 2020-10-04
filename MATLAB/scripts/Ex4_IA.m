%% IA vs LQR

close all
[m_sv, l_sv] = size(G);
u = u0(end-l_sv+1:end);
r = r_vec(end - m_sv+1:end);
[u_inf1, e_inf1, y_inf1, impr1, iteration_number1, error_history1] = RIA(G,d, b,r, u, 1);
xmax = iteration_number1;
if(xmax>0)
    xlim([0, xmax])
end
legend('$e_k$', 'interpreter', 'latex');

%% Plot results LQR vs IA 
close all
hold on
y_LQR = reshape(y_sv, [length(y_sv)/m, m]);
ref = reshape(r, [length(r)/m, m]);
y = reshape(y_inf1{end}, [length(y_inf1{end})/m, m]);

plot(0:length(y_LQR(:,1))-1, y_LQR(:,1));
plot(0:length(ref(:,1))-1, ref(:,1), 'LineWidth', 1.8);
plot(0:length(y(:,1))-1, y(:,1));
legend('LQR', 'Reference', 'IA');
%xlim([0, N-1])
hold off
