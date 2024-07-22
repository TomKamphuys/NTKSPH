function Ynm = spherical_harmonic(n, m, az, el)

  a = sqrt((2*n+1)./(4*pi) .* factorial(n-m) ./ factorial(n+m));

  % azimuthal changing part of spherical harmonics
  e = exp(1j*m.*az);

  % elevation dependend part of spherical harmonics
  l = legendre(n, cos(el))';
  l = l(:,abs(m)+1);
  if m < 0
      % legendre function for negative m from [2], eq. (6.31)
      l = l * (-1).^-m .* factorial(n+m)./factorial(n-m);
  end

  % get spherical harmonics
  Ynm = a .* l .* e;

endfunction
