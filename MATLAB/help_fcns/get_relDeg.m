function [rel_deg] = get_relDeg(A,B,C,D)
    
    rel_deg = 0;
    [m,l] = size(D);
    cont = 0;
    if all(all(D==0))
        cont = 1;
        rel_deg = rel_deg + 1;
    end
    
     A_mult = eye(length(A));
    while cont 
         M = C*A_mult*B;
         if all(all(M == 0))
             A_mult = A_mult*A; 
             rel_deg = rel_deg + 1;
         else
             cont = 0;
         end
    end

end

