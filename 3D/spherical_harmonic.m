function Y = spherical_harmonic(n, m, phi, theta)

  % TODO MPOT het lijkt dat de schaling nog niet goed is en het teken voor oneven n ook niet

  isoddm = mod(m,2) == 1;
  isnegm = m < 0;

  m = abs(m);

  P = legendre(n, cos(theta), 'norm');
  C = 1/sqrt(2*pi);
  P = P(abs(m)+1,:)';

  if isnegm
      E = sin(m*phi);
  else
      E = cos(m*phi);
  end

  Y = C * P .* E;

  if isoddm
      Y = -Y;
  end

endfunction
