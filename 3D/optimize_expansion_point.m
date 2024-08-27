function optimize_expansion_point()

  options = optimset();

  ind = 500;

%  for ind = 1000:100:1500
    disp(ind);

    corr = [0, 0, 0];
    error_dB = test_center_fit(corr, ind);
    disp(error_dB);

    x = fminsearch(@test_center_fit, corr, options, ind);
    disp(x);

    error_dB = test_center_fit(x, ind);
    disp(error_dB);

%    X(ind, :) = x;

%  endfor

endfunction
