%% 
% Returns system with non-zero D term for uniform reative degree corresp.
% 4.4.1
function [Anew,Bnew,Cnew,Dnew, Nnew] = get_non0D_system(A,B,C,D, N)
    rel_deg = get_relDeg(A,B,C,D);
    [m,l] = size(D);
   
    
    Anew = A;
    Bnew = B;
    if rel_deg > 0
        Dnew = C*A^(rel_deg - 1)*B;
        Cnew = C*A^rel_deg; 
        Nnew = N - rel_deg;
    else
        Cnew = C;
        Dnew = D;
        Nnew = N;
    end
    % u0new = u0(1:l*Nnew,1)


end

