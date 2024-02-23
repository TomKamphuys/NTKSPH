function [x, y, z] = conc_half_cyl(z_range, z_pnts, th_pnts, radii)
% CONC_HALF_CYL
%  Return measurement grid with points on concentric half cylindrical
%   surfaces centered to the z-axis. Grid points have non-negative
%   x-coordinates.
%
if nargin == 0
    z_range = [-0.50 0.50];
    z_pnts = 25;
    th_pnts = 25;
    radii = [0.25 0.30];
elseif nargin ~= 4
    error('Invalid number of input arguments.');
end
[z_single, th] = meshgrid(linspace(z_range(1), z_range(end), z_pnts), ...
                          linspace(-0.5*pi, 0.5*pi, th_pnts));
               
n_radii = numel(radii);
n_angs = numel(th);
x = zeros(n_radii*n_angs, 1);
y = zeros(n_radii*n_angs, 1);
z = zeros(n_radii*n_angs, 1);

for indx = 1:n_radii
    x(1+(indx-1)*n_angs:indx*n_angs) = radii(indx)*cos(th(:));
    y(1+(indx-1)*n_angs:indx*n_angs) = radii(indx)*sin(th(:));
    z(1+(indx-1)*n_angs:indx*n_angs) = z_single(:);
end

end