function yn = spherical_yn(nu, z)
% Compute the spherical Bessel function of the second kind y_nu(z)
%
% If both 'nu' and 'z' are length n and m vectors, a matrix will be
%  returned:
%  [ yn(nu_1, z_1)  yn(nu_2, z_1)  yn(nu_3, z_1) ... yn(nu_n, z_1) ]
%  [ yn(nu_1, z_2)  yn(nu_2, z_2)  yn(nu_3, z_2) ... yn(nu_n, z_2) ]
%  [      ...            ...            ...      ...      ...      ]
%  [ yn(nu_1, z_m)  yn(nu_2, z_m)  yn(nu_3, z_m) ... yn(nu_n, z_m) ]
%
cols = size(nu, 2);
rows = size(z, 2);
z_m = repmat(z, cols, 1)';
nu_m = repmat(nu, rows, 1);
yn = sqrt(pi ./(2 * z_m)) .* bessely(nu_m + 0.5, z_m);