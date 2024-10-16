function Xcir = ffcir(X)
[r,c,n] = size(X);
% Create the DFT matrix F_n and identity matrix I_m
F_n = dft_matrix(n);
I_r = eye(r);
I_c = eye(c);

% Compute Kronecker products
F_n_I = kron(F_n, I_r);
F_c_I = kron(F_n, I_c);
Xcir = F_n_I * bcir(X) * F_c_I';