function [spl, phase] = take_virtual_measurement(CD, phi, theta, r, freqs)

  N = calcNmax(size(CD, 1));

  temp = 273.15 + 20;
  kr = calc_kr(r, freqs, temp);

  angular_part = calc_angular_part(phi, theta, N);

  for ind = 1:length(kr)
    sph_hn1 = calc_radial_part(kr(ind), N) .* angular_part;
    p_recon = sph_hn1 * CD(:, ind);
    spl(:, ind) = dB_SPL(p_recon);
    phase(:,ind) = angle(p_recon)/pi*180;
  endfor

endfunction
