%return G and d in supervector form respective given system (A, B, C, D),
%zero initial condition x0 and number of iterations N,
%such that y = Gu + d in supervector form

function [G,d] = get_G(A, B, C, D, x0, N)
G = sparse(kron(eye(N+1),D));
d = eye(length(A))*x0;
tmp_mult_A = eye(length(A));%containes the last A^i, such that A^(i+1) can be calculated as A*tmp_mult_A
for i = 1:N
    G = G + sparse(kron(diag(ones(1,N-i+1),-i),C*tmp_mult_A*B));  
    
    tmp_mult_A = A*tmp_mult_A;
    d = [d; tmp_mult_A*x0];
end
d = kron(eye(N+1),C)*d;
end

