function hn2 = spherical_hn2(nu, z)
% Compute the spherical Hankel function of the second kind h^2_nu(z)
%
% If both 'nu' and 'z' are length n and m vectors, a matrix will be
%  returned:
%  [ h^2(nu_1, z_1)  h^2(nu_2, z_1)  h^2(nu_3, z_1) ... h^2(nu_n, z_1) ]
%  [ h^2(nu_1, z_2)  h^2(nu_2, z_2)  h^2(nu_3, z_2) ... h^2(nu_n, z_2) ]
%  [      ...             ...             ...       ...      ...       ]
%  [ h^2(nu_1, z_m)  h^2(nu_2, z_m)  h^2(nu_3, z_m) ... h^2(nu_n, z_m) ]
%
cols = size(nu, 2);
rows = size(z, 2);
z_m = repmat(z, cols, 1)';
nu_m = repmat(nu, rows, 1);
hn2 = sqrt(pi ./(2 * z_m)) .* besselh(nu_m + 0.5, 2, z_m);