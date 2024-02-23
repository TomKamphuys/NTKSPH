function [x, y, z, freqs, p_meas] = read_measurements(data_dir)
% READ_MEASUREMENTS
%  Read measurement data from files in sub-directory 'data_dir':
%   'coordinates.csv'
%   'frequencies.csv'
%   'measurements.csv'
%
% If 'data_dir' is not given or is an empty string, read from current
%  directory.
if nargin == 0
    data_dir = strcat('.', filesep);
elseif isempty(strtrim(data_dir))
    data_dir = strcat('.', filesep);
else
    if ~exist(data_dir, 'dir')
        error(strcat('Data directory "', data_dir, '" not found.'));
    end
end

coords = csvread(strcat(data_dir, filesep, 'coordinates.csv'));
x = coords(:, 1)';
y = coords(:, 2)';
z = coords(:, 3)';
freqs = csvread(strcat(data_dir, filesep, 'frequencies.csv'));
p_meas = csvread(strcat(data_dir, filesep, 'measurements.csv'));

end