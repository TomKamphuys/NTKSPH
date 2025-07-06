function [spl, angle, freqs] = horizontal_directivity(CD, freqs)

  distance = 10; % m

  n = 100;
  angle = linspace(-pi, pi, n);
  x_recon = distance .* cos(angle);
  y_recon = distance .* sin(angle);
  z_recon = zeros(size(y_recon));

  [phi, theta, r] = cart2sph_phys(x_recon', y_recon', z_recon');

  [spl, ~] = take_virtual_measurement(CD, phi, theta, distance, freqs);

endfunction


