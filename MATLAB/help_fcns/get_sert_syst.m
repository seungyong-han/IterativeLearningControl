function [A,B,C,D, Delta] =  get_sert_syst(L)
[Lsys, Delta] = lftdata(L);

dim = length(Lsys) - length(Delta);%Must be square

[M,N] = size(Lsys);

A = Lsys(1:dim, 1:dim);
B = Lsys(1:dim, dim+1:M);
C = Lsys(dim+1:M, 1:dim);
D = Lsys(dim+1:M, dim+1:N);

end