function CD = spherical_fit(r, theta, z, p, freq, Nmax)

  [x, y, z] = cyl2cart(r, theta, z);

  % Temperature is 20 deg C.
  temp = 273.15 + 20;

  % To make things simpler, simulate 1 frequency only
  freqs = freq;
  omega = 2*pi*freqs;

  % Simulate the measurements
  p_meas = p(:,round(freqs/10));
  x_meas = x;
  y_meas = y;
  z_meas = z;

  [phi_meas, theta_meas, r_meas] = cart2sph(x_meas, y_meas, z_meas);
  theta_meas = pi/2 - theta_meas;

  [PSI_mat, Nmax] = sph_PSI_mix(r_meas, theta_meas, phi_meas, omega, Nmax, temp);

  CD_vec = PSI_mat \ p_meas;

  CD = CD_vec(1:2:end);

endfunction


