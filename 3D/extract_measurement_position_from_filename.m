function [r, theta, z] = extract_measurement_position_from_filename(filename)
	
  position = strrep(filename, '(', '');
  position = strrep(position, ')', '');
  position = strsplit(position, ',');

  r = str2num(position{1})/1000; % convert to meters
  theta = deg2rad(str2num(position{2}));
  z = str2num(position{3})/1000; % convert to meters

endfunction