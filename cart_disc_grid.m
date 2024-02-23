function [x, y, z] = cart_disc_grid(plane, n_r, n_ang, r_span, ang_span)
% CART_DISC_GRID
%  Create a flat disc shaped plotting grid for surface plotting.
%  The grid may be on the xy, yz, or xz plane.
%  Return three (3) 2D arrays of x, y, z mesh coordinates defining the
%  disc surface and is intended for 3D mesh surface plots.
%  Inputs:
%   plane - valid choices are 'xz' (default), 'xy', or 'yz'
%   n_r - number of points in the radial direction
%   n_ang - number of points in the angular direction
%   r_span - row vector as [inner_radius  outer radius]
%   ang_span - row vector as [start_angle  end_angle] in radians
%     If plane is 'xz' or 'yz', ang_span specifies the span of polar
%     angles with the north pole being 0.
%
if nargin == 0
    plane = 'xz';
    n_r = 21;
    n_ang = 91;
    r_span = [1 2];
    ang_span = [0.0, 2*pi];
elseif nargin == 1
    n_r = 21;
    n_ang = 91;
    r_span = [1 2];
    ang_span = [0.0, 2*pi];
elseif nargin ~= 5
    error('Invalid number of input parameters.');
end

[m_r, m_ang] = meshgrid(linspace(r_span(1), r_span(end), n_r), ...
                        linspace(ang_span(1), ang_span(end), n_ang));

if strcmpi(plane, 'xz') || strcmpi(plane, 'zx')
   [x, y, z] = sph2cart(0, pi/2-m_ang, m_r);
elseif strcmpi(plane, 'yz') || strcmpi(plane, 'zy')
   [x, y, z] = sph2cart(pi/2, pi/2-m_ang, m_r);
elseif strcmpi(plane, 'xy') || strcmpi(plane, 'yx')
   [x, y, z] = sph2cart(m_ang, 0, m_r);
else
    error('Invalid value for parameter %s',plane);
end

end

