function axlim = plot_extent(x, y, z)
% PLOT_EXTENT
%  Determine the suitable plot limits given the coordinates 'x, y, z'.
%  'axlim' is the max absolute value of 'x', 'y', and 'z', round up
%   to the next 1 decimal place.
%
axlim = max([max(reshape(abs(x), 1, [])) max(reshape(abs(y), 1, [])) ...
             max(reshape(abs(z), 1, []))]);
axlim = 0.1 * (floor(10 * axlim + 1e-6) + 1) ;
end

