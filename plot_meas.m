function plot_meas(x, y, z, freqs, p_meas, f_bin, ax)
% PLOT_MEAS
%  Plot the measurements at frequency bin 'f_bin' to axis 'ax'.
%
if ~ishandle(ax)
    error('%s is not a valid handle of a graphics object.',n_plots);
end

f = freqs(f_bin);
str_f = sprintf('Bin # %d  f = %.0f Hz', f_bin, f);
scatter3(ax, x, y, z, abs(p_meas(:, f_bin)), abs(p_meas(:, f_bin)));
ax.Title.String = str_f;

if exist('params.mat', 'file')
    var_info = who('-file', 'params.mat');
    if ismember('cmap', var_info)
        load('params.mat', 'cmap');
        colormap(ax, cmap);
    end
    if ismember('caxlims', var_info)
        load('params.mat', 'caxlims');
        caxis(caxlims);
    end
    
end

axis equal;
axlim = plot_extent(x, y, z);
ax.XLim = [-axlim axlim];
ax.YLim = [-axlim axlim];
ax.ZLim = [-axlim axlim];
campos([1, 1, 1]);

end