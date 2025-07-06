function [x, y, z] = sph2cart_phys(phi, theta, r)

  % phi is the angle in the xy plane (phi = 0 when along x axis)
  % theta is the angle containing z (theta = 0 when along z axis)
  % phi and theta are different in the iso standard, which we don't use, but matlab/octave does
  phi_iso = pi/2 - theta;
  theta_iso = phi;
  [x, y, z] = sph2cart(theta_iso, phi_iso, r);

endfunction
