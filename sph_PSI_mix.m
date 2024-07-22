function [PSI_mat, Nout] = sph_PSI_mix(r, theta, phi, omega, N, temp)
% SPH_PSI_MIX
%  Compute the mixed field spherical wave expansion function matrix PSI.
%  The function 'AKsh' from AKtools is used to calculate the spherical
%   harmonics.
%  Inputs: 'r', 'theta', and 'phi' are the coordinates.
%          'omega' is frequency in rad/s.
%          'N' is the highest order of the coefficients to compute. When
%           not given, it is calculated from the number of measurement
%           locations (i.e. length(theta)) so that the resulting matrix
%           PSI is just over-determined (i.e. more rows than columns)
%          'temp' is temperature in K, default 293.15 K (= 20 degC)
%  Output: Computed spherical wave expansion function matrix PSI in
%           complex numbers. The size of the matrix is: M rows by
%           2*(N+1)**2 columns, where M = # of coordinates.
%


R_air = 287.058;

c = sqrt(1.4 * R_air * temp);
kr = r * omega / c;
n_rows = size(theta, 2);
n_cols = (N + 1)^2;
PSI_mat = zeros(n_rows, n_cols);

for n = 0:N
    hn1 = spherical_hn1(n, kr);
    jn = spherical_jn(n, kr);
    for m = -n:n
        j = 2*(n^2 + n + m) + 1;
        sph_harm = spherical_harmonic(n, m, phi, theta);
        PSI_mat(:, j) = sph_harm .* hn1;
        PSI_mat(:, j+1) = sph_harm .* jn;
    end
end

Nout = N;

end
