function PSI_mat = sph_PSI_mix(r, theta, phi, freqs, N, temp)
  % SPH_PSI_MIX
  %  Inputs: 'r', 'theta', and 'phi' are the coordinates.
  %          'freqs' is frequency in 1/s.
  %          'N' is the highest order of the coefficients to compute. When
  %           not given, it is calculated from the number of measurement
  %           locations (i.e. length(theta)) so that the resulting matrix
  %           PSI is just over-determined (i.e. more rows than columns)
  %          'temp' is temperature in K, default 293.15 K (= 20 degC)
  %  Output: Computed spherical wave expansion function matrix PSI in
  %           complex numbers. The size of the matrix is: M rows by
  %           2*(N+1)**2 columns, where M = # of coordinates.
  %


  kr = calc_kr(r, freqs, temp);

  radial_out = calc_radial_part(kr, N);
  radial_in = calc_radial_part_in(kr, N);
  angular = calc_angular_part(phi, theta, N);

  PSI_mat_in = angular .* radial_in;
  PSI_mat_out = angular .* radial_out;

  PSI_mat = [PSI_mat_out PSI_mat_in];

endfunction
