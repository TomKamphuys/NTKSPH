function [p, r, theta, z, f] = read_nfs_measurements(subdir)
  % read_nfs_measurements reads nfs measurements (plain wave files) containing
  % impulse responses measured at certain positions. The position in encoded in
  % the filename. It expects the data to be in a subdir of the Measurements directory.
  % It uses MATAA for some processing.


  unit = 'V';
  directory = ['Measurements/' subdir '/'];
  files = glob([directory '/*.wav']);

  for ind = 1:numel(files)
    [~, file] = fileparts(files{ind});

    [hh, fs] = audioread([directory file '.wav']);
    t = (1:length(hh))/fs;
    bla(ind) = sum(abs(hh(1:40)));
    [mag, phase, f, unit_mag] = mataa_IR_to_FR(hh,fs,[],unit);

    pt = to_pressure(mag) .* exp(i.*deg2rad(phase));

    p(ind, :) = pt';

    [r(ind), theta(ind), z(ind)] = extract_measurement_position_from_filename(file);

  endfor

% This code could be used in case only half a sphere is measured and you wan to mirror the measurements
%  p = [flipud(p); p(2:end,:)];
%  angles = deg2rad([flipud(-angles); angles(2:end)]);

endfunction
