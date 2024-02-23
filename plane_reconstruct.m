function [X, Y, out_ref, out_recon] = plane_reconstruct(x_size, y_size, amount, pt_srcs, CD_vec, omega, Nmax, temp)

  x = linspace(-x_size, x_size, amount);
  y = linspace(-y_size, y_size, amount);

  [X, Y] = meshgrid(x, y);

  x_recon = X(:)';
  y_recon = Y(:)';
  z_recon = zeros(size(y_recon));

  % Convert the reconstruction points coordinates from Cartesian to polar
[phi_recon, theta_recon, r_recon] = cart2sph (x_recon, y_recon, z_recon);
theta_recon = pi/2 - theta_recon;

PSI_recon = sph_PSI_mix(r_recon, theta_recon, phi_recon, omega, Nmax, temp);

p_recon = PSI_recon * CD_vec;

out_recon = dB_SPL(reshape(p_recon, size(X)));

out_ref = reshape(dB_SPL(sim_meas_cart(x_recon, y_recon, z_recon, omega, pt_srcs)), size(X));


endfunction
