function [out_recon, angle, freqs] = horizontal_directivity(CD, freqs)

  distance = 10; % m
  temp = 273.15 + 20;
  R_air = 287.058;

  N = calcNmax(size(CD, 1));

  n = 100;
  angle = linspace(-pi, pi, n);
  x_recon = distance .* cos(angle);
  y_recon = distance .* sin(angle);
  z_recon = zeros(size(y_recon));

  [phi, theta, r] = cart2sph_phys(x_recon', y_recon', z_recon');

  kr = calc_kr(r, freqs, temp);

  angular_part = calc_angular_part(phi, theta, N);
  for ind = 1:length(kr)
    sph_hn1 = calc_radial_part(kr(ind), N) .* angular_part;
    p_recon = sph_hn1 * CD(:, ind);
    out_recon(:, ind) = dB_SPL(p_recon);
  endfor

endfunction


