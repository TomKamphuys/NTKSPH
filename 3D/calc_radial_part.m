function out = calc_radial_part(kr, N)

  out = zeros(length(kr), (N + 1)^2);
  for n = 0:N
    hn1 = spherical_hn1(n, kr);
    for m = -n:n
      j = n^2 + n + m + 1;
      out(:, j) = hn1;
    endfor
  endfor
endfunction
