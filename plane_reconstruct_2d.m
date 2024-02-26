function [X, Y, out_ref, out_recon] = plane_reconstruct_2d(x_size, y_size, amount, pt_srcs, CD_vec, omega, Nmax, temp)

  x = linspace(-x_size, x_size, amount);
  y = linspace(-y_size, y_size, amount);

  [X, Y] = meshgrid(x, y);

  x_recon = X(:)';
  y_recon = Y(:)';

  % Convert the reconstruction points coordinates from Cartesian to polar
  [theta, r] = cart2pol(x_recon, y_recon);

  PSI_recon = sph_PSI_mix_2d(r, theta, omega, Nmax, temp);

  p_recon = PSI_recon * CD_vec;

  out_recon = dB_SPL(reshape(p_recon, size(X)));

  out_ref = reshape(dB_SPL(sim_meas_cart_2d(x_recon, y_recon, omega, pt_srcs)), size(X));


endfunction
