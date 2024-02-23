%% ----- Step 1 -----
%{
  Initial release: 2021-04-18
  Updated:         2021-04-22

  Read measurement data.
  'data_dir': Subdirectory with the measurement data files. Current
              directory if not specified or empty string
%}
[x, y, z, freqs, p_meas] = read_measurements(data_dir);
fprintf('Number of measurement points: %d\n', length(x))
fprintf('Number of frequencies       : %d\n', length(freqs))
% % Plot some of the measurement data
% if ~isempty(meas_plots)
%     fig1 = grp_plot_meas(x, y, z, freqs, p_meas, meas_plots);
%     drawnow;
% end