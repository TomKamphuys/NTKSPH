function [THETA, Z, out_ref, out_recon] = cylwall_reconstruct(y_size, amount, CD_full, omega, Nmax, temp)

  theta = linspace(-pi, pi, amount);
  z = linspace(-y_size, y_size, amount);

  [THETA, Z] = meshgrid(theta, z);

  [X, Y, Z] = cyl2cart(3, THETA, Z);

  x_recon = X(:)';
  y_recon = Y(:)';
  z_recon = Z(:)';

  % Convert the reconstruction points coordinates from Cartesian to polar
  [phi_recon, theta_recon, r_recon] = cart2sph_phys(x_recon, y_recon, z_recon);

  PSI_recon = sph_PSI_mix(r_recon', theta_recon', phi_recon', omega, Nmax, temp);

  p_ref = PSI_recon * CD_full;

  n = size(CD_full, 1);
  p_recon = PSI_recon(:, 1:n/2) * get_outgoing_coefficients(CD_full);  % SFS reconstruction!

  out_ref   = dB_SPL(reshape(p_ref, size(X)));
  out_recon = dB_SPL(reshape(p_recon, size(X)));

endfunction
