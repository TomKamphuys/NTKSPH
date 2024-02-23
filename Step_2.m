%% ----- Step 2 -----
%{
  Initial release: 2021-04-18
  Updated:         2021-05-04

  The fitting process to solve for the spherical wave expansion function 
   coefficients works at one frequency at a time.
%}

%apply offset correction based on step 1b

% Convert Cartesian coordinates 'x', 'y', 'z' to spherical coordinates
%  'r', 'theta', 'phi'

%===APPLY OFFSET CORRECTION================
% [phi, theta, r] = cart2sph (x, y, z);
[phi, theta, r] = cart2sph(x-offsetposition(1,1), y-offsetposition(1,2), z-offsetposition(1,3));

% MATLAB returns 'theta' as elevation angles. Need to convert to polar
%  angles.
theta = pi/2 - theta;
% Optionally to convert the range of 'phi' from [-pi pi] to [0 2pi]
% phi = phi + 2*pi*(phi < 0);

% With 'Nmax' = -1, 'sph_PSI_mix()' will determine the order Nmax so that
%  the PSI matrix is not under-determined. The PSI matrices for the lower
%  orders are subsets of the PSI matrix calculated to 'Nmax', and we
%  don't need to calculate them separately.
% If 'Nmax' isn't already defined, set to -1.
if ~exist('Nmax', 'var')
    Nmax = -1;
end

% For this example, use the frequency bin 'f_bin' to reconstruct
omega = 2*pi * freqs(f_bin);

% Time how long it takes to calculate the PSI matrix.
fprintf('Compute measurement grid PSI matrix to order %d ...\n', Nmax)
tic;
[PSI_mat, Nmax] = sph_PSI_mix(r, theta, phi, omega, Nmax, temp);
fprintf('Measurement grid PSI matrix: ')
fprintf('Nmax = %d, compute time = %.3f, ', Nmax, toc)
fprintf('size %d X %d\n', size(PSI_mat))