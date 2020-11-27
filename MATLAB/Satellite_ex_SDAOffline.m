G_star = sparse((R\eye(length(R)))*G'*Q);

beta = beta*2/max(abs(eig(G*G_star)));  
K0 = G_star; 

%define initial parameter
u = u0;
e0 = r - G*u - d; 
e = e0;
cont = 1;
iteration_number = 0;
error_history = [norm(e0)];

betaGK0 = beta*G*K0;%calculate once
% sequencies setup
u_inf = [];
y_inf = [];
cell_nb = 0;
while cont
    iteration_number = iteration_number + 1;
    
    %calculate u_new and e_new with update rules
    u_new = u + beta*K0*e;
    u_new(u_new<-10) = -10;
    u_new(u_new>10) = 10; 
    %add saturation
    
    
    e_new = r - G*u_new - d;%(eye(length(e0)) - betaGK0)*e;
    
    
    if(any(isnan(u_new)) || any(isnan(e_new)))
        disp('NAN')
    end
    
    
    
    %Termination criterion
    if norm(e - e_new)<10^-6
        cont = 0;
    end
   
    %print current imprivement
    if(mod(iteration_number, 10000) == 0)
        disp("currErrordiff:")
        norm(e_new - e)
    end
    
    u_inf = [u_inf, u_new]; 
    y_inf = [y_inf, G*u_new + d]; 

    
    e = e_new; 
    u = u_new;
    error_history = [error_history, norm(e)];

end



impr = norm(e0)/error_history(end);
