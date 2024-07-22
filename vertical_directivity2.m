function [out_recon, angles, freqs] = vertical_directivity2(CD, freqs)

  distance = 10; % m
  % Temperature is 20 deg C.
  temp = 273.15 + 20;
  Nmax = calcNmax(size(CD, 1));

  n = 100;
  angles = linspace(-pi, pi, n);
  x_recon = distance .* cos(angles);
  y_recon = zeros(size(x_recon));
  z_recon = distance .* sin(angles);

  % Convert the reconstruction points coordinates from Cartesian to polar
  [phi_recon, theta_recon, r_recon] = cart2sph(x_recon, y_recon, z_recon);
  theta_recon = pi/2 - theta_recon;

  R_air = 287.058;
  c = sqrt(1.4 * R_air * temp);
  factor = 2*pi/c

  for ind = 1:length(freqs)

    kr = r_recon * factor * freqs(ind);

    PSI_recon = sph_PSI_out_fast(theta_recon, phi_recon, kr, Nmax);
    p_recon = PSI_recon * CD(:, ind);
    out_recon(:, ind) = dB_SPL(p_recon);

  endfor

endfunction

function PSI_mat = sph_PSI_out_fast(theta, phi, kr, N)

  hn1 = first_hankel_fast(kr, N);
  sph_harm = spherical_harmonics_fast(theta, phi, kr, N);

  PSI_mat = sph_harm .* hn1;

endfunction

function hn1 = first_hankel_fast(kr, N)
    for n = 0:N
      hn1(:, n) = spherical_hn1(n, kr);
    endfor

    % TODO nog een repmat maar hoeveel precies???

endfunction

function sph_harm = spherical_harmonics_fast(theta, phi, kr, N)

  n_rows = size(theta, 2);
  n_cols = (N + 1)^2;
  PSI_mat = zeros(n_rows, n_cols);

  for n = 0:N
      for m = -n:n
          j = (n^2 + n + m) + 1;
          sph_harm(:, j) = spherical_harmonic(n, m, phi, theta);
      end
  end

endfunction
