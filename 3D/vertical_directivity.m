function [spl, angle, freqs] = vertical_directivity(CD, freqs)

  distance = 10; % m
  R_air = 287.058;
  temp = 273.15 + 20;
  N = calcNmax(size(CD, 1));

  n = 100;
  angle = linspace(-pi, pi, n);
  x_recon = distance .* cos(angle);
  y_recon = zeros(size(x_recon));
  z_recon = distance .* sin(angle);

  [phi, theta, r] = cart2sph_phys(x_recon', y_recon', z_recon');

  [spl, ~] = take_virtual_measurement(CD, phi, theta, distance, freqs);

endfunction

