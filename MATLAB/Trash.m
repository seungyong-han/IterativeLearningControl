
if length(G(1,:)) == length(G(:, 1)) && abs(det(G))~=0 
        [u_inf, e_inf, y_inf, impr,iteration_number, error_history]  = RIA(G,d, beta,r, u0, R, Q, do_plot);
        disp('Use right inverse model algorithm')
        return 
end


%% Satellite Expl -- doesnt work :( 

clear all
clc

J1 = 1;
J2 = .1;
k = .091; 
b = .0036; 

A = [0 1 0 0; -k/J2 -b/J2 k/J2 b/J2; 0 0 0 1; k/J1 b/J1 -k/J1 -b/J1];
B1 = [0 1 0 0]';
B2 = [0 0 0 1/J1]';

C = [1 0 0 0]; 
D11 = 0;
D12 = 0;

% %With disturbance: 
% B = [B1 B2]; 
% D = [D11 D12];

% Without disturbance
B = B2;
D = D12;


G = ss(A, B, C, D);


%[K,CL,gamma] = hinfsyn(G,1,1)

s = zpk('s');
Kn = 7.9212*(s + 0.1818)*(s^2 - 0.2244*s + 0.8981);
Kd = (s^2 + 3.899*s + 4.745)*(s^2 + 1.039*s + 3.395);
K = Kn/Kd;

Lc = (1+G*K)^-1;
isstable(Lc)
G = c2d(G, 1e-12);
K = c2d(K, 1e-12);

Ld = (1+G*K)^-1;
isstable(Ld)

[NumK, DenK] = tfdata(K)
NumK = NumK{1}
DenK = DenK{1}
[NumG, DenG] = tfdata(G)
NumG = NumG{1}
DenG = DenG{1}


%Find discrete controller for G 
b = .0396;
s = zpk('s');
G = (.036*(s + 25.28))/(s^2*(s^2 + b*s + 1)) + 1;
G = c2d(G, 1e-5)



hinfsyn(G, 1 ,1 ,1)
K = hinfsyn(G, 1 ,1 ,1)
tf(K)
tfdata(K)
[NumK, DenK] = tfdata(K)
NumK = NumK{1}
DenK = DenK{1}
[NumG, DenG] = tfdata(G)
NumG = NumG{1}
DenG = DenG{1}


L = (1+G*K)^-1
isstable(L)

    


%weights
s = zpk('s');
w1 = (.5*s + .433)/(s+.00433);
w2 = .1;

P1dof = [-1 G; 0 1; -1 G];

Kn = 7.9212*(s + 0.1818)*(s^2 - 0.2244*s + 0.8981);
Kd = (s^2 + 3.899*s + 4.745)*(s^2 + 1.039*s + 3.395);
K = Kn/Kd;
w1 = (.5*s + .433)/(s+.00433);
w2 = .1;
w3 = 1; 
P1dof = [w1 0 0; 0 w2 0; 0 0 1]*P1dof;
%P1dof = c2d(P1dof, 1);




%%

clc 
clear all 

s = zpk('s')
b = .0396;
G = (.036*(s + 25.28))/(s^2*(s^2 + b*s + 1));
G = c2d(G, 1e-5);

[A, B, C, D] = ssdata(G);

poles = [.3, -.4, .3+.3*i, .3-.3*i]
F = place(A, B, poles)
poles = [.3, -.4, .3+.3*i, .3-.3*i]
J = place(A',C',poles)'

Ak = A + B*F + J*C+ J*D*F;
Bk = J;
Ck = F;
[m ,l] = size(D);
Dk = zeros(m,l)

K = ss(Ak, Bk, Ck, Dk, 1e-5)

L = (1 +G*K)^-1;
isstable(L)





%%
%calculate transform matrix
    M = full(G*K0);
    [T,eig_M] = eig(M);
    T = T/norm(T^-1)^2;
    
    %!!!Works only for G*K0 diagonalizable!!!
    
    
    
    
    u = u0;
    e0 = r - G*u - d; 
    e = T^-1*e0;
    cont = 1;
    iteration_number = 0;
    error_history = [norm(e0)];
    input_history = [norm(u0)];
    
   
    dim_M = length(M);
    
    beta_1 = ones(10, 1)*1/norm(eig_M(1:10, 1:10));
    beta_2 = ones(dim_M-20,1)*1/norm(eig_M(11:dim_M-10, 11:dim_M-10));
    beta_3 = ones(10,1)*1/norm(eig_M(dim_M-9:dim_M, dim_M-9:dim_M));
    beta_ = beta;
    
    beta = T*diag([beta_1; beta_2; beta_3]);
%     alpha_ = 1/norm(eig_M*beta*T);

%%
     e = T^-1*e0;
     u = u0;
     %betaK0 = beta*K0;
     %betaGK0 = beta*G*K0;
     alpha_ = 1/norm(max(eig(eig_M*beta*T)));
     while cont
        iteration_number = iteration_number + 1;
        u_new = u + K0*beta*T*e;;
        e_new = (eye(length(e0)) - alpha_*eig_M*beta_*T)*e;
        

        if norm(e - e_new)<10^-10
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

    u_inf = u_new;
    e_inf = norm(e_new);
    y_inf = G*u_new + d;
    impr = norm(e0)/e_inf;

    if do_plot
        plot(0:iteration_number, error_history);
    end
    
%%

%% Test SDA (Alg. 7.3, 7.4, 7.5)
%   GG* singular
N = 10;
A = diag([-.2, -.2]);
B = [1 0 3; 0 2 1]; 
C = [.1 .5; -.4 .1; .3 -.8];
D = diag([1,4,0]);
x0 = [0, 0]';


%l-input m-output system; in supervector description l*(N+1) - input
%m*(N+1)-output
l = length(B(1,:));
m = length(C(:,1));

%simulate system response
sys = ss(A,B,C,D,1);

[G, d] = get_G(A, B, C, D, x0, N);
u0 = repmat(0.1,l*(N+1),1);
%r = G*G'*5*rand(12,1) + G*u0 + d;%to ensure that e0 is in image(GG*)
r = 5*rand(m*(N+1),1);
%r = [0 0 -1, zeros(1,9)]';%Element of ker(GG*) %TODO: Fixme Dimensions
%must agree

R = eye(l*(N+1));
Q = eye(m*(N+1));
beta = .2

do_plot = 1;
[u_inf1, e_inf1, y_inf1, impr1, iteration_number1] = SDA(G,d, beta,r, u0, R, Q,do_plot)
figure 
[u_inf2, e_inf2, y_inf2, impr2, iteration_number2] = SDA_suppression_of_evs(G,d, r, u0, R, Q,do_plot)
figure
[u_inf3, e_inf3, y_inf3, impr3, iteration_number3] = SDA_beta_as_evs(G,d, r, u0, R, Q,do_plot)



%% Uncertainty 

u0 = 5*rand(l*(N+1),1);
r = 5*rand(m*(N+1),1);
e0 = r - G.Nominal*u0;%FIXME d in get_G 
hold on
for t = 1:10
    evec = [];
    L0 = usample(L); %Wi bestimme ich aus welchem Range darf a sein wrt rausgefundener Norm fuer Delta? 
    e_new = e0;
    e = 2*e0;
    if max(abs(eig(L0)))<1
    while(norm(e_new - e)>1e-6)
        e = e_new;
        e_new = L0*e;
        evec = [evec; norm(e)];
    end
    else 
        disp('||L0||>=1')
    end
    
    plot(length(evec), evec)
    
end

hold off


Epsilon weighted K0
eps = .5;
G_eps = get_G(eps^2*A, eps^2*B, C, D, x0, N);

K0_eps = G_eps';
beta = .5/norm(G_eps.Nominal)^2;
L = (eye(m*(N+1)) - beta*G*K0_eps);
[Lsys, Delta] = lftdata(L);
blockStr = [-12 0; 4 4]
[bounds] = mussv(Lsys, blockStr)


LMIs Versuche
u0 = repmat(0.1,l*(N+1),1);
r = 5*rand(m*(N+1),1);

R = eye(l*(N+1));
Q = eye(m*(N+1));

beta = .1;
do_plot = 1;
do_print = 1;


%For IA 
K0 = inv(G);
L = (eye(m*(N+1)) + G*K0)



L = ss(L, zeros(length(L),1),zeros(1, length(L)),0,ts)

[Lsys, Delta] = lftdata(L);
[A_L, B_L, C_L, D_L] = ssdata(Lsys);

n_L = length(A_L);
l_L = length(B_L(1,:));
m_L = length(C_L(:,1));%m_L = l_L



setlmis([])
X = lmivar(1,[n_L,1]);

%For zero convergence 
% A1 = [eye(n_L) zeros(n_L, l_L); A_L B_L];
% A2 = [C_L D_L; zeros(l_L, n_L) eye(l_L)];
% A1I = A1'*[-eye(n_L) zeros(n_L); zeros(n_L) eye(n_L)]*A1;
% 
% lmiterm([1 1 1 0], -eye(n_L) + A_L'*A_L);
% lmiterm([1 1 1 Q], C_L', C_L);
% lmiterm([1,1,2,0], A_L'*B_L);
% lmiterm([1,1,2,Q], C_L', D_L);
% lmiterm([1,2,1,0], B_L'*A_L);
% lmiterm([1,2,1,Q], D_L', C_L);
% lmiterm([1,2,2,0], B_L'*B_L);
% lmiterm([1,2,2,Q], D_L', D_L);
% lmiterm([1,2,2,Q], -1,1);

%for non-zero convergence 

lmiterm([1 1 1 X], -1, 1)
lmiterm([1 1 1 X], A_L', A_L)
lmiterm([1 2 1 X], B_L', A_L)
lmiterm([1 1 2 X], A_L', B_L)
lmiterm([1 2 2 X], B_L', B_L)


lmis = getlmis; 
[tmin,xfeas] = feasp(lmis);



NON-Zero-Robust-Convergence
X = sdpvar(length(L))
cnstr = [X >= 0];
Ls = ssdata(usample(L))
cnstr = cnstr + [Ls'*X*Ls - Ls <= 0]
solvesdp(cnstr);    

[Lsert, Delta] = lftdata(L)
get(Delta)


[A_Lsert, B_Lsert, C_Lsert, D_Lsert] = ssdata(Lsert);


a0 = .1;
a = sdpvar(length(D_Lsert), 1); 
Delta = diag(a);
Delta = uncertain(Delta)
X = sdpvar(length(L));
Y = sdpvar(length(D_Lsert));
% A_Lsert = Lsert(1,1);
% B_Lsert = Lsert(1,2:13);
% C_Lsert = Lsert(2:13, 1);
% D_Lser = Lsert(2:13, 2:13);

Lusert = A_Lsert + B_Lsert*Delta*C_Lsert;
cnstr = [X >= 0];
cnstr = [(eye(13) + D_Lsert*Delta)*Y == eye(13)]
cnstr = cnstr + [Lusert'*X*Lusert - Lusert <= 0]
cnstr = cnstr + [Delta <= .4*eye(13)]
cnstr = cnstr + [Delta >= -.2*eye(13)]
solvesdp(cnstr)


