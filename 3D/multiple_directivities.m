function multiple_directivities(CD, freqs, steps)

  distance = 10; % m

  n = 101;
  angle = linspace(-pi, pi, n);
  x_recon = distance .* cos(angle);
  y_recon = distance .* sin(angle);
  z_recon = zeros(size(y_recon));


  hoeken = linspace(0, pi/2, steps);

  for ind = 1:steps
    rotation_angle = hoeken(ind);
    x = x_recon;
    y = y_recon*cos(rotation_angle) - z_recon*sin(rotation_angle);
    z = y_recon*sin(rotation_angle) + z_recon*cos(rotation_angle);

    [phi, theta, r] = cart2sph_phys(x', y', z');

    [spl, ~] = take_virtual_measurement(CD, phi, theta, distance, freqs);

    spl = spl - repmat(spl(ceil(n/2),:), [n, 1]);

    figure(ind)
    pcolor(freqs, angle, spl)
    shading flat
    set(gca,'xscale','log');
    colormap jet
    colorbar
    caxis([-40 0]) % adjust for actual limits
    title(sprintf('Rotation angle %f degrees', rotation_angle/pi*180))
  endfor


endfunction
