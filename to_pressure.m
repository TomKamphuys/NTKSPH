function out = to_pressure(in)
%  dB_out = 20 * log10(abs(x) / 20e-6);

  out = 10.^(in ./ 20) .* 20e-6;


endfunction
