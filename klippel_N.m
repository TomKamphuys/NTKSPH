function out = klippel_N(f)

  freq = [0 500 600 700 800 900 1000 1200 1400 1600 1800];
  N = [5 5 6 7 8 9 10 11 12 13 14];

  out = interp1(freq, N, f, 'previous');

endfunction
