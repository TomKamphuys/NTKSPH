function plot_2d_functions(N, k)
% plot 2d functions

angle = linspace(-pi, pi, k);
value = cos(angle);

figure
hold on

if N == 0;
  lengte = 1;
else
  lengte = ceil(sqrt(N+1));
end


for n = 0:N

    l = myLegendre(n, value);

    subplot(lengte, lengte, n+1)
%    polar(angle, abs(l));
    color_polar_plot(angle, l);
    title(sprintf('n = %d', n));
    grid off
endfor

axis equal;

endfunction
