function [Anew, Bnew, Cnew, Dnew] = get_stable_ss(G, P, K)
[Ak, Bk, Ck, Dk] = ssdata(K); 
[Acl,Bcl,Ccl,Dcl] = ssdata(lft(P,K)); 
[A,B,C,D] = ssdata(G); 
l = length(B(1,:)); 
m = length(C(:, 1)); 
n = length(A); 

Anew = Acl; 

Cnew = (eye(length(D)) - D*Dk)^-1*[C D*Ck];
Dnew = (eye(length(D) - D*Dk))^-1*D; 

Bnew = blkdiag(B, Bk)*...
    ([eye(length(Dk(:,1)), length(D(1,:))),  - Dk;...
     - D, eye(length(D(:,1)), length(Dk(1,:)))])^-1*...
    [zeros(m,l); D] + [B; zeros(length(Ak),l)];
end

