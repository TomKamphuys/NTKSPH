function [p, r, theta, z, f] = read_nfs_measurements()

  unit = 'V';
  decimation_amount = 1;
  date = '21062024-spherical'
  directory = ['Measurements/' date '/'];
  files = glob([directory '/*.wav']);

  for ind = 1:numel(files)
    [~, file] = fileparts(files{ind});

    position = strrep(file, '(', '');
    position = strrep(position, ')', '');
    position = strsplit(position, ',');

    [h, fs] = audioread([directory file '.wav']);
    t = (1:length(h))/fs;
%    [t_start,t_rise] = mataa_guess_IR_start(h,fs);
%    [hh,th] = mataa_signal_crop(h,fs,t_start-t_rise,t(end));
    hh = h;
    t = (1:length(h))/fs;
    [mag,phase,f,unit_mag] = mataa_IR_to_FR(hh,fs,[],unit);

    pt = to_pressure(mag) .* exp(i.*deg2rad(phase));
    f = f(1:decimation_amount:end);
    pt = pt(1:decimation_amount:end);

    p(ind,:) = pt';

    r(ind) = str2num(position{1})/1000;
    theta(ind) = str2num(position{2})/180*pi;
    z(ind) = str2num(position{3})/1000;

  endfor

%  p = [flipud(p); p(2:end,:)];
%  angles = deg2rad([flipud(-angles); angles(2:end)]);

endfunction
