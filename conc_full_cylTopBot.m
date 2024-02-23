function [x, y, z] = conc_full_cylTopBot(z_range, z_pnts, th_pnts, radii, r_pnts_cap)
% CONC_FULL_CYL
%  Return measurement grid with points on concentric full cylindrical
%   surfaces centered to the z-axis.
%

if nargin == 0
    z_range = [-0.50 0.50];
    z_pnts = 15;
    th_pnts = 72;
    radii = [0.25 0.30];
    r_pnts_cap = 5;
elseif nargin ~= 5
    error('Invalid number of input arguments.');
end

% %%========== comment for function===========
% clear all;

% z_range = [-0.75 0.75];
% z_pnts = 1;
% th_pnts = 4;
% radii = [1.0 2 ];
% r_pnts_cap=1;

% %%========== comment for function===========


[z_single, th] = meshgrid(linspace(z_range(1), z_range(end), z_pnts), linspace(-pi, pi-((2*pi)/th_pnts), th_pnts));
               
n_radii = numel(radii);
n_angs = numel(th);
n_pnts_cap = th_pnts*r_pnts_cap;

x = zeros(n_radii*n_angs + (n_pnts_cap), 1);
y = zeros(n_radii*n_angs + (n_pnts_cap), 1);
z = zeros(n_radii*n_angs + (n_pnts_cap), 1);

for indx = 1:n_radii
    x(1+(indx-1)*n_angs:indx*n_angs) = radii(indx)*cos(th(:));
    y(1+(indx-1)*n_angs:indx*n_angs) = radii(indx)*sin(th(:));
    z(1+(indx-1)*n_angs:indx*n_angs) = z_single(:);
end

% %add points for top and bottom caps
capradii = linspace( (radii(1)/(r_pnts_cap+1)) , radii(1)-(radii(1)/(r_pnts_cap+1)),r_pnts_cap);

count=n_angs*n_radii+1;
zi(1)=z_range(1);
zi(2)=z_range(end);
radiioffset(1)=0;
radiioffset(2)=abs(radii(1)-radii(end));

for radiii = 1:n_radii
    for zii= 1:2
        for capradiii = 1:r_pnts_cap
            for thi= 1:th_pnts  

   x(count)= capradii(capradiii)*cos(th(thi));
   y(count)= capradii(capradiii)*sin(th(thi));
       if zii==1
           z(count)= zi(1)-radiioffset(radiii); 
       else
           z(count)= zi(2)+radiioffset(radiii) ;
       end; 

   
   count=count+1;
             end
         end
     end
end

% end
