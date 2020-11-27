%% Apply ILC algorithms

N = 0;

%l-input m-output system; in supervector description l*(N+1) - input
%m*(N+1)-output
l = length(B(1,:));
m = length(C(:,1));

beta = .6;
cont = 50; 
Nmax = 3e4; 

%% IA
u_inf = {}; 


Gstr=  {}; 
dstr = {}; 

%Apply the algorithm for count first elements 
%and add Gstr, dstr
[An,Bn,Cn,Dn, Nnew] = get_non0D_system(A - B*F, B, C-D*F, D, cont);
diff = cont - Nnew; 
for i = 1:cont
    i = i - diff;
    u = u_sv(1:i+1); 
    r = r_vec(1:i+1); 
    
    [G, d] = get_G(An, Bn, Cn, Dn, x0, i); 
      
    [u_inf, e_inf, y_inf, impr,iteration_number, error_history] = RIA(G,d, beta,r, u, 0);%SDA(G,d, beta,r, u, 1, 1, 0);
    u_sv(1:i+1) = u_inf{end}; 
    i = i + diff; 
end


N = cont - diff; 
while N<Nmax
    N = N+1;
    
    u = u_sv(N+1 - cont:N);%SISO system 
    r = r_vec(N+1 - cont:N); 
    x0 = A*x0 + B*u_sv(N+1-cont); 
    [G, d] = get_G(An, Bn, Cn, Dn, x0, Nnew); 

    
    [u_inf, e_inf, y_inf, impr,iteration_number, error_history] = RIA(G,d, beta,r, u, 0);
    if(all(r>0))
    close all
    plot(1:length(r), r);
    hold on 
    plot(1:length(y_inf{end}), y_inf{end}, '--');
    hold off
    end
%     if(N==2499)
%         disp('here')
%     end
    
    u_sv(N+1 - cont:N) = u_inf{end}; 
    if mod(N, 1000) == 0
        disp(['Iteration ', num2str(N), ' of ',num2str(Nmax)]); 
    end
end
%%
x = zeros(length(A),1); 
for i = 1:Nmax
    x = A*x + B*u_sv(i); 
    
    y(i) = C*x + D*u_sv(i); 
    
end

plot(1:length(y), y); 
hold on
plot(1:length(r_vec), r_vec)
hold off


