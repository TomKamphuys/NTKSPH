function [phi, theta, r] = cart2sph_phys(x, y, z)

  [phi, theta, r] = cart2sph(x, y, z);
  theta = pi/2 - theta;

endfunction
