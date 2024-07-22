function db = evaluate(r, theta, phi, omega, CD, temp)

  N = calcNmax(size(CD, 1));

  for ind = 1:size(r)
    PSI = sph_PSI_out(r(ind), theta(ind), phi(ind), omega, N, temp);
    p = PSI * CD;

    db(ind) = dB_SPL(p);

  endfor

endfunction
