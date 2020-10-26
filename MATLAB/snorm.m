function nrm = snorm(M)
nrm = sqrt(max(abs(eig(M*M'))));
end

