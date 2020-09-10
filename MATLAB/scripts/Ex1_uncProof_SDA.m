%% Robustness proof for SDA
a = ureal('a', 4, 'PlusMinus', [-3, 3]);
A = [2 1; a 3];
A = A - B*F;
ts = 1e-3;

usys = ss(A, B, C, D, ts);
[A, B1, B2, C1, D11, D12, C2, D21, D22, Delta] = get_expanded_system(usys);
[G1, d1] = get_G(A, B1, C1, D11, x0, N-1);
G2 = get_G(A, B2, C1, D12, x0, N-1);
[G3,d2] = get_G(A, B1, C2, D21, x0, N-1);
G4 = get_G(A, B2, C2, D22, x0, N-1);

%K = .1*G^-1; 
K = 1/norm(G)^2*G';

AL = eye(length(G1)) + G1;
BL = G2;
CL = G3*K;
DL = G4;

X = sdpvar(length(AL));
%X = [2 0; 0 .5];

%X = eye(2);
X = kron([-1,0; 0 1], X);

P = sdpvar(length(DL));
Q = sdpvar(length(DL));
P = kron([1 0; 0 -1], P) + kron([0 1; 1 0], Q);


M1 = [eye(length(AL)), zeros(length(AL), length(BL(1,:))); AL, BL]
M2 = [CL, DL; zeros(length(DL),length(CL(1,:))), eye(length(DL))];
cnstr = [M1'*X*M1  + M2'*P*M2<=0]
cnstr = cnstr + [Q' + Q == 0]


%cnstr = cnstr + [-X + A'*X*A <= 0, A;*X*B]

sol = solvesdp(cnstr) 

