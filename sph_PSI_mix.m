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
if nargin < 4
    error('Too few input parameters.');
elseif nargin == 4
    N = -1;
    temp = 293.15;
elseif nargin == 5
    temp = 293.15;
elseif nargin > 6
    error('Too many input parameters.');
end

if N < 0
    N = floor(sqrt(length(theta)/2) - 1);
% elseif N > floor(sqrt(length(theta)/2) - 1)
%     fprintf('N = %d will result in underdetermined system.\n', N)
%     fprintf('Length(theta) = %d\n', length(theta))
%     N = floor(sqrt(length(theta)/2) - 1);
%     fprintf('Reduced N to %d\n', N)
end

if exist('params.mat', 'file')
    var_info = who('-file', 'params.mat');
    if ismember('R_air', var_info)
        load('params.mat', 'R_air');
    else
        R_air = 287.058;    % Specific gas constant for air
    end
else
    R_air = 287.058;
end

c = sqrt(1.4 * R_air * temp);
kr = r * omega / c;
n_rows = size(theta, 2);
n_cols = 2 * (N + 1)^2;
PSI_mat = zeros(n_rows, n_cols);

for n = 0:N
    hn1 = spherical_hn1(n, kr);
    jn = spherical_jn(n, kr);
    for m = -n:n
        j = 2 * (n^2 + n + m) + 1;
        sph_harm = AKsh(n, m, rad2deg(phi), rad2deg(theta));
        PSI_mat(:, j) = sph_harm .* hn1;
        PSI_mat(:, j+1) = sph_harm .* jn;
    end
end

if nargout > 1
    Nout = N;
end

end