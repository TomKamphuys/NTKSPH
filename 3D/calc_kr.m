function kr = calc_kr(r, freqs, temp)

  R_air = 287.058;
  c = sqrt(1.4 * R_air * temp);
  omega = 2*pi*freqs;
  kr = r .* omega ./ c;

endfunction
