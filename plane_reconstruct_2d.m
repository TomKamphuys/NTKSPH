function [X, Z, out_ref, out_recon] = plane_reconstruct_2d(x_size, z_size, amount, pt_srcs, CD_vec, omega, Nmax, temp)

  x = linspace(-x_size, x_size, amount);
  z = linspace(-z_size, z_size, amount);

  [X, Z] = meshgrid(x, z);

  x_recon = X(:)';
  y_recon = zeros(size(x_recon)); %Y(:)';
  z_recon = Z(:)';

  % Convert the reconstruction points coordinates from Cartesian to polar
  [phi_recon, theta_recon, r_recon] = cart2sph (x_recon, y_recon, z_recon);
  theta_recon = pi/2 - theta_recon;

  PSI_recon = sph_PSI_mix_2d(r_recon, theta_recon, phi_recon, omega, Nmax, temp);

  p_recon = PSI_recon * CD_vec;

  out_recon = dB_SPL(reshape(p_recon, size(X)));

  out_ref = reshape(dB_SPL(sim_meas_cart(x_recon, y_recon, z_recon, omega, pt_srcs)), size(X));


endfunction
