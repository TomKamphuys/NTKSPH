function fig1 = grp_plot_meas(x, y, z, freqs, p_meas, n_plots)
% GRP_PLOT_MEAS
%  Plot a group of the measurement data.
%   'x', 'y', 'z' coordinates for all measurement points. These
%     coordinates are the same for each of the individual plots.
%   'freqs' is the list of the measurement frequencies.
%   'p_meas' is a 2D array of the measured pressure amplitudes.
%   'n_plots' is the list of frequency bins to plot.
%
if length(n_plots) < 1
    error('Array length %s of must be at least 1',n_plots)
end

if length(n_plots) > length(freqs)
    n_plots = 1:length(freqs);
end

n_rows = floor(sqrt(length(n_plots)));
if mod(length(n_plots), n_rows) > 0
    n_cols = floor(length(n_plots) / n_rows) + 1;
else
    n_cols = length(n_plots) / n_rows;
end

aspect = n_cols / n_rows;
fig1 = figure('Name', 'Selected Measurement Data', ...
              'Position', [100 190 850*aspect 700]);
icount = 1;
for n = n_plots
    ax = subplot(n_rows, n_cols, icount);
    plot_meas(x, y, z, freqs, p_meas, n, ax);
    icount = icount + 1;
end

end