%Angeblich falsch

function [u_inf, e_inf, y_inf, impr, iteration_number, error_history] = FD(G,d, Kc,beta,r, u0, do_plot)



S = (eye(length(G))  + G*Kc)\eye(length(G));
T = 1 - S;
beta = beta*2/norm(T)/norm(T');

l = length(G(1,:));%input dimension (length of the row)
m = length(G(:,1));%output dimension (length of the column)

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

M1 = beta*Kc*S*T';
M2 = beta*(T*T');

while cont
    u_new = u + M1*e;
    e_new = (eye(m) - M2)*e;
    
    
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
end



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