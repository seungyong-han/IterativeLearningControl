%% Apply ILC algorithms

iteration_number = 0;
N = 0;
u_inf = [];
e_inf = [];
y_inf = [];
impr = [];
act_trValue = [];
tic
colNb = 1;
plots = {};
plots{1} = [];
plots{2} = [];
plots{3} = [];
plots{4} = [];
plots{5} = [];
count = 100;
Astr = {eye(length(A))}; 
for i = 1:count
    [G,d] = get_redG(A-B*F, B, C-D*F, D, x0, i);
    Gstr{i} = G;
    dstr{i} = d; 
    Astr{i+1} = Astr{i}*A; 
end

Gmdl = Gstr;
dmdl = dstr; 



%l-input m-output system; in supervector description l*(N+1) - input
%m*(N+1)-output
l = length(B(1,:));
m = length(C(:,1));

beta = .3;

do_plot = 0;
do_print = 0;
u_hist = {};

N = count; 

r = r_vec(1:N+1); 
u0 = u_sv(1:N+1); 
[u_inf, e_inf, y_inf, impr,iteration_number, error_history] = RIA(Gstr{end},dstr{end}, beta,r, u0, 0);

u_sv(1:N+1) = u_inf{end}; 


while N < Nmax - 1
    N = N + 1; 
    u0 = u_sv(N - count+1:N+1); 
    r = r_vec(N - count+1:N+1);
    [u_inf, e_inf, y_inf, impr,iteration_number, error_history] = RIA(Gstr{end},dstr{end}, beta,r, u0, 0);

end




while N < Nmax-1
    
    
    N = N + count;

    %need to calculate new x0 
%     for i = 1:count
%         [G,d] = get_redG(A-B*F, B, C-D*F, D, x0, i);
%         Gstr{i} = G;
%         dstr{i} = d; 
%     end
%     Gmdl = Gstr;
%     dmdl = dstr; 

    
    u0 = u(N - count + 1:N+1, :);
    r = r_vec(N - count + 1:N+1, :);
    for i = 1:count
        u_curr = u0(1:i+1);    
        r_curr = r(1:i+1); 
        G = Gstr{i};
        d = dstr{i};
        e = r_curr - G*u_curr -d; 
        
        error_history = norm(e);
        R = eye(l*(i+1));
        Q = eye(m*(i+1));
        cont = 1;
        iteration_number = 0;
        
        K0 = pinv(G); 
        while cont
            iteration_number = iteration_number + 1;
            u_new = u_curr + beta*K0*e;
        
            e_new = r_curr - Gmdl{i}*u_new - dmdl{i}; %for demonstration non-stability

            if norm(e - e_new)<10^-6
                cont = 0;
            end
            error_history = [error_history, norm(e_new)];
            if(mod(iteration_number, 100) == 0)
                disp("currErrordiff:")
                norm(e_new - e)
            end
            u_curr = u_new; 
            e = e_new;
            
            if i == count
                if iteration_number == 20
                    u_hist{1} = [u_hist{1}, u_curr]; 
                end
                if iteration_number == 40
                    u_hist{2} = [u_hist{2}, u_curr]; 
                end
                if iteration_number == 80
                    u_hist{3} = [u_hist{3}, u_curr]; 
                end
            end
        end
        u0(1:i+1) = u_curr; 
    end
    
    if mod(N, 1000) == 0
        disp(['Iteration ', num2str(N), ' of ', num2str(Nmax)]); 
    end

end 
ell_time = toc;

