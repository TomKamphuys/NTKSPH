function view_modes(N)
  nr = 100;

  [theta, phi] = meshgrid(linspace(0, pi, nr), linspace(-pi, pi, nr));

  theta = theta(:);
  phi = phi(:);

  figure
  hold on
  for n = 0:N
    for m = -n:n
      value = spherical_harmonic(n, m, phi, theta);
      teken = sign(value);
      r = abs(value);

      [x, y, z] = sph2cart_phys(phi, theta, r);

      surf(reshape(x, [nr, nr]), reshape(y, [nr, nr])+m, reshape(z, [nr, nr])-2*n, reshape(teken, [nr, nr]))

    end
  end

  axis equal
  shading flat
  colormap summer
  colorbar
  camlight
  axis off

endfunction
