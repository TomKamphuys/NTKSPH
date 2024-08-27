function [error_dB, error_percentage] = calc_error(p_meas, PSI, CD_vec)

  p_mod = PSI * CD_vec;

  error_percentage = sum((abs(p_mod - p_meas)).^2) / sum((abs(p_meas)).^2) * 100;
  error_dB = 10*log10(error_percentage/100);

endfunction
