function [x, mse] = lstsq_solve(A, b)
% LSTSQ_SOLVE
%  Use MATLAB backslash to compute the least squares fit, where:
%   x = b \ A    solves for 'x' in the equation    A x = b
%  The mean square error of the least squares residuals 'mse' is
%   returned along with the vector 'x'
%
% Turn off the ill-conditioned matrix warning
% dA = decomposition(A, 'CheckCondition', false);
%
% x = dA \ b;

x = A \ b;

mse = norm(A*x - b)^2 / length(b);
end
