function balloon_plot(CD, freqs, index)

  temp = 273.15 + 20;
  R_air = 287.058;
  nr = 100;
  distance = 3;

  N = calcNmax(size(CD, 1));

  [theta, phi] = meshgrid(linspace(0, pi, nr), linspace(-pi, pi, nr));

  kr = calc_kr(distance, freqs, temp);

  angular_part = calc_angular_part(phi(:), theta(:), N);
  sph_hn1 = calc_radial_part(kr(index), N) .* angular_part;
  p = sph_hn1 * CD(:, index);
  out = dB_SPL(p);

  r = reshape(abs(p), [nr, nr]);

  [x, y, z] = sph2cart_phys(phi, theta, r);

  surf(x, y, z, r)
  axis vis3d
  shading interp

endfunction
