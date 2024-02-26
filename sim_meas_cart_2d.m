function p = sim_meas_cart_2d(x, y, omega, pt_srcs)
% SIM_MEAS_CART
%  Simulate the complex pressure amplitudes at locations 'x', 'y'
%  at frequency 'omega' due to point sources in 'pt_srcs'.
%  'pt_srcs' is a cell array in Cartesian coordinates in the form:
%    {[amplitude1, x1, y1, z1, phase1] [...] ... [...]}
%  'temp' is temperature in K.
% Outputs: Array of complex pressure amplitudes, same dimension as 'x'
%
if ismatrix(x)
    x_is_matrix = 1;
    x_size = size(x);
    x = x(:);
    y = y(:);
end

if exist('params.mat', 'file')
    var_info = who('-file', 'params.mat');
    if ismember('temp', var_info)
        load('params.mat', 'temp');
    else
        temp = 293.15;
    end
    if ismember('R_air', var_info)
        load('params.mat', 'R_air');
    else
        R_air = 287.058;
    end
else
    temp = 293.15;
    R_air = 287.058;
end
c = sqrt(1.4 * temp * R_air);   % Speed of sound

p = complex(zeros(size(x)));
for n = 1:length(pt_srcs)
    ctr_x = pt_srcs{n}(2);
    ctr_y = pt_srcs{n}(3);
    r = sqrt((x - ctr_x).^2 + (y - ctr_y).^2);
    % Limit r to be at least 1e-5 to avoid singularity
    r = max(r, repmat(1e-5, size(r)));
    p = p + (pt_srcs{n}(1)./r) .* exp(1i*(r*omega/c + pt_srcs{n}(4)));
end

if x_is_matrix
    p = reshape(p, x_size);
end

end
