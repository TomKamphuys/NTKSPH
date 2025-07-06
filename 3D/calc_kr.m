function kr = calc_kr(r, freqs, temp)
  % fix: r en freqs kunnen beide vectoren zijn. matrix van maken

  R_air = 287.058;
  c = sqrt(1.4 * R_air * temp);
  omega = 2*pi*freqs;
  kr = r .* omega ./ c;

endfunction
