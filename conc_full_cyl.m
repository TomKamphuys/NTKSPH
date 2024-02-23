function [x, y, z] = conc_full_cyl(z_range, z_pnts, th_pnts, radii)
% CONC_FULL_CYL
%  Return measurement grid with points on concentric full cylindrical
%   surfaces centered to the z-axis.
%
if nargin == 0
    z_range = [-0.50 0.50];
    z_pnts = 25;
    th_pnts = 25;
    radii = [0.25 0.30];
elseif nargin ~= 4
    error('Invalid number of input arguments.');
end
[z_single, th] = meshgrid(linspace(z_range(1), z_range(end), z_pnts), linspace(-pi, pi-((2*pi)/th_pnts), th_pnts));
               
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