%% ----- Step 0 -----
%{
  Initial release: 2021-04-18
  Updated:         2021-05-04

  Set some of the simulation parameters.
  Also creates the file 'params.mat' to pass some of the parameters to
   functions called later in the sesson.
%}
temperature=20;
temp = 273.15 + temperature;   % Temperature in K
R_air = 287.058;   % Specific gas constant for air
c = sqrt(1.4 * R_air * temp);%soundspeed at temperature in m/s
data_dir = '';   % Subdirector with the measurement data
meas_plots = [1 2];  % Plot the measurements of these bins.
                           % No plotting if array is empty.
f_bin = 2;  % Frequency bin used for the reconstruction
Nmax = 12;  % Maximum expansion order. '-1' will maximize the maximum
            %  expsnion order while keeping the PSI matrix from being
            %  under-determined based on the number of measurements.
UseOfffsetCorrection =1; %set if calculated offset correction is used

% Save these parameters to 'params.mat' for functions that may not have
%  access to the variables in the main workspace.
cmap = 'jet';  % Colormap used for 3D surface and scatter plots
caxlims = [-40 6];  % Colormap range limits
save('params.mat', 'cmap', 'caxlims', 'temp', 'R_air','c');