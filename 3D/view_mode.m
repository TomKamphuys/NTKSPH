function view_mode(n, m)
  nr = 100;

  [theta, phi] = meshgrid(linspace(0, pi, nr), linspace(-pi, pi, nr));

  el = theta(:);
  az = phi(:);

  value = harmonicY(n, m, el, az, 'type', 'real');
  teken = sign(real(value));
  r = abs(value);
%  r = real(spherical_harmonic(n, m, az, el));

  r = reshape(r, [nr, nr]);

  [x, y, z] = sph2cart_phys(phi, theta, r);

  surf(x, y, z, reshape(teken, [nr, nr]) .* r)
  axis equal
  shading flat
  colormap jet

endfunction
