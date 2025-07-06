function write_measurements(hor, ver)

  formatSpec = '%f %f %f\n';


  for ind = 1:length(hor.angle)
    filename = sprintf('hor %d.txt', hor.angle(ind));
    values = [hor.f(ind,:)' hor.spl(ind,:)' hor.phase(ind,:)']';
    fid = fopen(filename, 'w');
    fprintf(fid, formatSpec, values);
    fclose(fid);
  end

  for ind = 1:length(ver.angle)
    filename = sprintf('ver %d.txt', ver.angle(ind));
    fid = fopen(filename, 'w');
    values = [ver.f(ind,:)' ver.spl(ind,:)' ver.phase(ind,:)']';
    fprintf(fid, formatSpec, values);
    fclose(fid);
  end

endfunction
