function [theta, phi] = equi_dist_angles(n_approx)
% EQUI_DIST_ANGLES
%  Generate the angular coordinates of a grid with equi-distance
%   spacings in the polar ('theta') and azimuthal ('phi')
%   directions.
% 'n_approx' is the approximate number of grid points to be
%   generated. (The algorithm used may not be able to generate the
%   exact specified number of points, but only to a close
%   approximate.)
%
% The regular placement algorithm by Markus Deserno of Max-Planck-
%  Institut fur Polymerforschung, Ackermannweg 10, 55128 Mainz, Germany
%  is used to generate the angular coordinates.
%   https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
%
if nargin == 0
    n_approx = 1000;
elseif nargin ~= 1
    error('This function accepts only one integer input argument.');
end

theta = [];
phi = [];
a = 4 * pi / n_approx;
d = sqrt(a);
M_th = round(pi / d);
d_theta = pi / M_th;
d_phi = a / d_theta;

for m = 1:fix(M_th)
    th = pi * (m - 0.5) / M_th;
    M_ph = round(2 * pi * sin(th) / d_phi);
    for n = 1:fix(M_ph)
        theta = [theta; th];
        phi = [phi; 2 * pi * (n - 1) / M_ph];
    end
end

end
