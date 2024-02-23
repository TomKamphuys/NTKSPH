function [x, y, z] = planar(y_range, y_pnts, z_range, z_pnts, x_list)
% PLANAR
%  Return measurement grid with points on rectangular planar grids
%   parallel to YZ plane.
%
if nargin == 0
    y_range = [-0.50 0.50];
    y_pnts = 25;
    z_range = [-0.50 0.50];
    z_pnts = 25;
    x_list = [0.25 0.30];
elseif nargin ~= 5
    error('Invalid number of input arguments.');
end

[z, y] = meshgrid(linspace(z_range(1), z_range(end), z_pnts), ...
                  linspace(y_range(1), y_range(end), y_pnts));

npts_1 = numel(z);
z = repmat(z(:), length(x_list), 1);
y = repmat(y(:), length(x_list), 1);
x = zeros(size(z));

for indx = 1:length(x_list)
    x(1+(indx-1)*npts_1:indx*npts_1) = x_list(indx);
end

end

