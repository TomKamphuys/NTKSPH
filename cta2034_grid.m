function [r, theta, phi] = cta2034_grid()
% CTA2034_GRID
%  Return the coordinates of the ANSI/CTA-2034-A spinorama 70 points grid
%   in the spherical coordinate system.
%  'r' - row vector for the radii
%  'theta' - row vector for the polar angles (North pole = 0,
%            South pole = pi)
%  'phi' - row vector for the azimuthal angles (positive x-axis = 0,
%          positive y-axis = pi/2, negative x-axis= pi,
%          negative y-axis = 3pi/2)
%
r = ones(1, 70) * 2;
theta = [0:10:350  ones(1, 34)*90] * pi/180;
phi = [zeros(1, 36)  10:10:170  190:10:350] * pi/180;

end

