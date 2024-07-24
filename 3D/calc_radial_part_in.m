function out = calc_radial_part_in(kr, N)

  out = zeros(length(kr), (N + 1)^2);
  for n = 0:N
    jn = spherical_jn(n, kr);
    for m = -n:n
      j = n^2 + n + m + 1;
      sph_hn1(:, j) = jn;
    endfor
  endfor

endfunction
