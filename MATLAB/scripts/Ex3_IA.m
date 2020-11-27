%% Apply ILC algorithms

N = 0;

%l-input m-output system; in supervector description l*(N+1) - input
%m*(N+1)-output
l = length(B(1,:));
m = length(C(:,1));

beta = .6;

%% IA
u_inf = {}; 


Gstr=  {}; 
dstr = {}; 

%Apply the algorithm for count first elements 
%and add Gstr, dstr
[An,Bn,Cn,Dn, Nnew] = get_non0D_system(A, B, C, D, cont);
diff = cont - Nnew; 

u_sv = reshape(u_sv, [(Nmax+1)*m, 1]);
r_vec = reshape(r_vec, [(Nmax+1)*m, 1]);
x_sv2 = x0; 
for i = 1:cont - diff
    u = u_sv(diff*m+1:(i+1+diff)*m); 
    r = r_vec(diff*m+1:(i+1+diff)*m); 
    
    [G, d] = get_G(An, Bn, Cn, Dn, x0, i); 
      
    [u_inf, e_inf, y_inf, impr,iteration_number, error_history] =SDA(G,d, beta,r, u, 1, 1, 0);% RIA(G,d, beta,r, u, 0);%
    u_sv(diff*m+1:(i+1+diff)*m) = u_inf{end}; 
    
    x0 = A*x0 + B*u_sv(diff*m+1:(diff*m + 1)*m);
    x_sv2 = [x_sv2; x0]; 
end


N = cont - diff; 
%N = cont; 
while N<Nmax
    N = N+1;
    
    u = u_sv((diff + N-cont)*m + 1:(N + diff + 1)*m);
    r = r_vec((diff + N-cont)*m + 1:(N + diff + 1)*m);
    x0 = An*x0 + Bn*u_sv((diff + N - cont)*m + 1:(N+1 - cont+diff)*m);
    [G, d] = get_G(An, Bn, Cn, Dn, x0, Nnew); 
    x_sv2 = [x_sv2; x0]; 
    
    [u_inf, e_inf, y_inf, impr,iteration_number, error_history] = SDA(G,d, beta,r, u, 1, 1, 0);%RIA(G,d, beta,r, u, 0);
%     if(all(r>0))
%     close all
%     plot(1:length(r), r);
%     hold on 
%     plot(1:length(y_inf{end}), y_inf{end}, '--');
%     hold off
%     end
   
    
    u_sv((diff + N-cont)*m + 1:(N + diff + 1)*m) = u_inf{end}; 
    if mod(N, 1000) == 0
        disp(['Iteration ', num2str(N), ' of ',num2str(Nmax)]); 
    end
end
%%
x = zeros(length(A),1); 
for i = 1:Nmax
    x = A*x + Bn*u_sv(i); 
    
    y(i) = Cn*x + Dn*u_sv(i); 
    
end

plot(1:length(y), y); 
hold on
plot(1:length(r_vec), r_vec)
hold off


