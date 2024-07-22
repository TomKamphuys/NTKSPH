function out = fun(f, CD, r_recon, theta_recon, phi_recon, Nmax, temp)

  omega = 2*pi*f;

  PSI_recon = sph_PSI_out(r_recon, theta_recon, phi_recon, omega, Nmax, temp);
  p_recon = PSI_recon * CD(:, ind);
  out = dB_SPL(p_recon);

endfunction
