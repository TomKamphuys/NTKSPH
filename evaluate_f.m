function db = evaluate_f(theta, phi, CD, hn1)

  N = calcNmax(size(CD, 1));


  sph_harm = doe_iets(N, phi, theta);

  for ind = 1:size(CD, 2)
    PSI = sph_harm .* hn1(ind, :);
    p = PSI * CD(:, ind);

    db(ind) = dB_SPL(p);
  endfor

endfunction

function sph_harm = doe_iets(N, phi, theta)

n_rows = size(theta, 2);
n_cols = (N + 1)^2;
sph_harm = zeros(n_rows, n_cols);

for n = 0:N
    for m = -n:n
        j = (n^2 + n + m) + 1;
        sph_harm(j) = spherical_harmonic(n, m, phi, theta);
    end
end

endfunction

