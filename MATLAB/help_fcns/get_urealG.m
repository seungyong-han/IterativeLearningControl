%returns G and d in supervector for given uncertain system (A, B, C, D),
%zero initial condition x0 and number of iterations N,
%such that y = Gu + d in supervector forms

function [G,d] = get_urealG(A, B, C, D, x0, N)
l = length(D(1,:));
m = length(D(:,1));

G = [D, zeros(m, N*l)];

A_mult = {eye(length(A))};

for i = 2:N+1
    A_mult{i} = A_mult{i-1}*A;
end

%calculate G
for i = 1:N
   Gtmp = [];
   for  j = i:-1:1
       Gtmp = [Gtmp,C*A_mult{j}*B];
   end
   Gtmp = [Gtmp D];
   
   Gtmp = [Gtmp, zeros(m, (N) - i)];
   
   G = [G; Gtmp];
end

%calculate d
d = C*x0;
for i = 2:N+1
    d = [d; C*A_mult{i}*x0];
end
end
