function [r, theta, phi] = cta2034_halfgrid()
% CTA2034_HALFGRID
%  Return the coordinates of the front half of the ANSI/CTA-2034-A
%   spinorama 70 points grid in the spherical coordinate system.
%  'r' - row vector for the radii
%  'theta' - row vector for the polar angles (North pole = 0,
%            South pole = pi)
%  'phi' - row vector for the azimuthal angles (positive x-axis = 0,
%          positive y-axis = pi/2, negative x-axis= pi,
%          negative y-axis = 3pi/2)
%
r = ones(1, 37) * 2;
theta = [0:10:180  ones(1, 18)*90] * pi/180;
phi = [zeros(1, 19)  10:10:90  270:10:350] * pi/180;

end

