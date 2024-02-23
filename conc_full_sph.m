function [x, y, z] = conc_full_sph(n_approx, radii)
% CONC_FULL_SPH
%  Return a measurement grid consists of concentric full spherical
%   surfaces centered to the origin. The grid points on each surface
%   are approximately equally spaced. Angular coordinates of the points
%   on each spherical surface are the same as all other surfaces.
%   The number of concentric spherical surfaces is determined by the
%   number of elements in array 'radii'.
%  The algorithm does not guarantee it will generate the requested
%   number of point, but only a close number to it.
%  'n_approx' - approximate number of points to be generated.
%  'radii' - an array of radii for the spherical surfaces.
%
if nargin == 0
    n_approx = 1250;
    radii = [0.45 0.50];
elseif nargin ~= 2
    error('Invalid number of input arguments.');
end

[theta, phi] = equi_dist_angles(round(n_approx/length(radii)));

n_angs = numel(theta);
n_radii = numel(radii);
x = zeros(n_radii * n_angs, 1);
y = zeros(n_radii * n_angs, 1);
z = zeros(n_radii * n_angs, 1);

for indx = 1:n_radii
    x(1+(indx-1)*n_angs:indx*n_angs) = radii(indx)*cos(phi).*sin(theta);
    y(1+(indx-1)*n_angs:indx*n_angs) = radii(indx)*sin(phi).*sin(theta);
    z(1+(indx-1)*n_angs:indx*n_angs) = radii(indx)*cos(theta);
end

end