function [Ak, Bk, Ck, Dk, x0k] = get_reduced_system(A, B, C, D, x0, Order, reducer)
    G = ss(A,B,C,D, 1); 

    %Tranform initial condition
    [K,g,T,Ti]=balreal(G);
    x0k = T*x0;

    % Perform balanced truncation on LTI system
    System = G; % Define System to reduce

    % Create option set for balred command
    Options = balredOptions();

    % Compute reduced order approximation
    ReducedSystem = balred(System,Order,Options);


    [Ak, Bk, Ck, Dk, stk] = ssdata(ReducedSystem);
    x0k = x0k(1:length(Ak));
end