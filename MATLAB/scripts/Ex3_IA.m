%% Apply ILC algorithms

iteration_number = 0;
Nmax = 3e4+1;%TODO: change to N
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
count = 50;
[G, d] = get_G(A,B, C, D, x0, count);
while N < Nmax-1l
    
    
   
    N = N + count;
    
     
    %l-input m-output system; in supervector description l*(N+1) - input
    %m*(N+1)-output
    l = length(B(1,:));
    m = length(C(:,1));
    u0 = u(N - count + 1:N+1, :);
    


    
    R = eye(l*(count+1));
    Q = eye(m*(count+1));
    
    beta = .1;
    r_curr = r(N - count + 1:N+1, :);
    do_plot = 0;
    do_print = 0;
    [u_inf1, e_inf1, y_inf1, impr1, iteration_number1] = SDA(G,d, beta,r_curr, u0, R, Q,do_plot);
    u_inf = [u_inf; u_inf1{length(u_inf1)}];
    y_inf = [y_inf; y_inf1{length(y_inf1)}];
    e_inf = [e_inf; e_inf1{length(e_inf1)}];
    act_trValue = [act_trValue; r];
    impr = [impr; impr1];
    iteration_number = iteration_number + iteration_number1;
    if mod(N, 1e5)==0
    fprintf('Iteration %d of %d\n', N+1, Nmax);
    end

       q = 5;
    for t = 1:q
            plots{t} = [plots{t}; y_inf1{t}];
    end
    
    
end 
ell_time = toc;

