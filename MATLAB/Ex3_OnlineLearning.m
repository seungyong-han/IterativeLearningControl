%Define all the matrices G will be needed
Gstr = {};
dstr = {};
count = 50;

for i = 1:count
    [G,d] = get_redG(A-B*F, B, C-D*F, D, x0, i);
    Gstr{i} = G;
    dstr{i} = d; 
end

Gmdl = Gstr;%stays for the real plant -- eg y = Gmdl{i}*u + d equals  
            %the real measurements data at time i; here: = Gstr, but can be
            %replaced by perturbated model 

%% IA
beta = .1;

beta_max = 2;
beta = beta*beta_max;

u = u0(1:2*l);

r0 = r_vec;
r = r0(1:2*m);

for i = 1:count
    G = Gstr{i};
    d = dstr{i}; 

    K0 = pinv(full(G));
    %K0 = G\eye(length(G));
    K0 = sparse(K0);

    e0 = r - Gmdl{i}*u - d; 
    e = e0;
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
    end
    if i<N
        u = [u; u0(i+1:(i+1)*l)];
        r = r0(1:(i+2)*m); 
    end
end


if do_plot
    plot(0:iteration_number, error_history);
end
