function jn = spherical_jn(nu, z)
% Compute the spherical Bessel function of the first kind j_nu(z)
%
% If both 'nu' and 'z' are length n and m vectors, a matrix will be
%  returned:
%  [ jn(nu_1, z_1)  jn(nu_2, z_1)  jn(nu_3, z_1) ... jn(nu_n, z_1) ]
%  [ jn(nu_1, z_2)  jn(nu_2, z_2)  jn(nu_3, z_2) ... jn(nu_n, z_2) ]
%  [      ...            ...            ...      ...      ...      ]
%  [ jn(nu_1, z_m)  jn(nu_2, z_m)  jn(nu_3, z_m) ... jn(nu_n, z_m) ]
%
cols = size(nu, 2);
rows = size(z, 2);
z_m = repmat(z, cols, 1)';
nu_m = repmat(nu, rows, 1);
jn = sqrt(pi ./(2 * z_m)) .* besselj(nu_m + 0.5, z_m);