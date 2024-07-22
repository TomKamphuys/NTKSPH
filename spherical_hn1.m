function hn1 = spherical_hn1(nu, z)

  hn1 = sqrt(pi ./(2 * z)) .* besselh(nu + 0.5, 1, z);

endfunction
