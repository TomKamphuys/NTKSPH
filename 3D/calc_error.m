function [error_dB, error_percentage] = calc_error(r, theta, phi, p_meas, CD_vec, omega, Nmax, temp)

  PSI = sph_PSI_mix(r, theta, phi, omega, Nmax, temp);

  p_mod = PSI * CD_vec;

  error_percentage = sum((abs(p_mod - p_meas)).^2) / sum((abs(p_meas)).^2) * 100;
  error_dB = 10*log10(error_percentage/100);

endfunction
