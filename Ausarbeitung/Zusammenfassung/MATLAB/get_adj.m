function [G_star] = get_adj(R, G, Q)
G_star = sparse((R\eye(length(R)))*G'*Q);
end

