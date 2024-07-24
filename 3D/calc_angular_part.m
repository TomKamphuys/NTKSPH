function out = calc_angular_part(phi, theta, N)
  out = zeros(length(phi), (N + 1)^2);
  for n = 0:N
    for m = -n:n
      j = n^2 + n + m + 1;
      out(:, j) = spherical_harmonic(n, m, phi, theta);
    endfor
  endfor
endfunction
