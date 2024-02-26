function color_polar_plot(a, r)

  teken = (r >= 0);

  changes = find(diff(teken))+1;

  if isempty(changes)
    if (teken(1) > 0)
      polar(a, abs(r), 'b');
    else
      polar(a, abs(r), 'r');
    end

    return;
  end


  previous = 1;
  for ind = 1:length(changes)
    if (teken(changes(ind)-1) > 0)
      polar(a(previous:changes(ind)), abs(r(previous:changes(ind))), 'b');
      hold on;
    else
      polar(a(previous:changes(ind)), abs(r(previous:changes(ind))), 'r');
      hold on;
    end

    previous = changes(ind);
  end

  if (teken(changes(ind)+1) > 0)
    polar(a(previous:end), abs(r(previous:end)), 'b');
    hold on;
  else

    polar(a(previous:end), abs(r(previous:end)), 'r');
    hold on;
  end

  hold off;

endfunction
