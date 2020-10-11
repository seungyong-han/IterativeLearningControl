%Define all the matrices G will be needed
Gstr = {};
dstr = {};
for i = 1:N
    [G,d] = get_redG(A-B*F, B, C-D*F, D, x0, N);
    Gstr{i} = G;
    dstr{i} = d; 
end

Gmdl = Gstr;%stays for the real plant -- eg y = Gmdl{i}*u + d equals  
            %the real measurements data at time i

%% IA
beta = .1;


for i = 1:N
    G = Gstr{i};
    d = dstr{i}; 
    beta_max = 2;
    beta = beta*beta_max;
    beta = beta*2;
    K0 = pinv(full(G));
    %K0 = G\eye(length(G));
    K0 = sparse(K0);

    u = u0(1:i);
    e0 = r - G*u - d; 
    e = e0(1:i);
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
        e_new = r - Gmdl{i}*u_new - d; %for demonstration non-stability

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

