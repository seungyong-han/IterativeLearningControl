function [A, B1, B2, C1, D11, D12, C2, D21, D22, Delta] = get_expanded_system(usys)
    [sys, Delta] = lftdata(usys);
    
    [A, B, C, D] = ssdata(sys);
    
    dim = size(Delta);
    
    B1 = B(:, dim(1)+1:end);
    B2 = B(:, 1:dim(2));
    
    C1 = C(dim(1)+1:end, :);
    C2 = C(1:dim(1), :);
    
    
    D22 = D(end-dim(1)+1:end, end-dim(2)+1:end);
    D11 = D(1:end-dim(1), 1:end - dim(2));
    D12 = D(1:end-dim(1), end-dim(2)+1:end);
    D21 = D(end-dim(1)+1:end, 1:end - dim(2));

end