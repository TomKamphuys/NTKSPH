function [THETA, Z, out_ref, out_recon] = cylwall_reconstruct(y_size, amount, CD_vec, omega, Nmax, temp)

  theta = linspace(-pi, pi, amount);
  z = linspace(-y_size, y_size, amount);

  [THETA, Z] = meshgrid(theta, z);

  [X, Y, Z] = cyl2cart(3, THETA, Z);

  x_recon = X(:)';
  y_recon = Y(:)';
  z_recon = Z(:)';

  % Convert the reconstruction points coordinates from Cartesian to polar
[phi_recon, theta_recon, r_recon] = cart2sph(x_recon, y_recon, z_recon);
theta_recon = pi/2 - theta_recon;


PSI_recon = sph_PSI_mix(r_recon, theta_recon, phi_recon, omega, Nmax, temp);

p_ref = PSI_recon * CD_vec;

%CD_vec = zeros(81,1);
%CD_vec(11) = .1;
p_recon = PSI_recon(:, 1:2:end) * CD_vec(1:2:end, :);  % SFS reconstruction!


out_ref   = dB_SPL(reshape(p_ref, size(X)));
out_recon = dB_SPL(reshape(p_recon, size(X)));

endfunction
