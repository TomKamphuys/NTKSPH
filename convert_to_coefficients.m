function CD = convert_to_coefficients()

  load 21062024_measurements.mat;

  N = 8;
  temp = 273.15 + 20;

  z = z;
  r = r + 0.03;

  [x, y, z] = cyl2cart(r, theta, z);

  x_meas = x - 0.13; % acoustic center correction
  y_meas = y;
  z_meas = z - 0.09; % acoustic center correction

  [phi, theta, r] = cart2sph_phys(x_meas', y_meas', z_meas');

  sph_harm = calc_angular_part(phi, theta, N);

  kr = calc_kr(r, f, temp);

  for ind = 1:length(f)
    waitbar(ind/length(f));

    p_meas = p(:, ind);

    sph_hn1 = calc_radial_part(kr(ind), N) .* sph_harm; % TODO !!!! zowel in as outgoing stuk doen

    CD_vec = sph_hn1 \ p_meas;

    CD(:, ind) = CD_vec(1:2:end); % TODO nieuwe conventie is niet in en out om de beurt, maar als 2 volleidge blokken achter elkaar.
  endfor

endfunction

