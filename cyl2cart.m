function [x, y, z] = cyl2cart(r, t, z)

  x = r.*cos(t);
  y = r.*sin(t);
  z = z;

endfunction
