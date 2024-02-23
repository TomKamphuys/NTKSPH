function [eX,eb] = algebraicGPSequations(satPos, tVals)
% Calculates the position of a receiver given the time delay of arrival to
% at least 4 satellites in 3D (or at least 3 satellites in 2D).
% INPUTS: satPos is Cartesian location of satellites (in 1D, 2D or 3D)
%         tVals is the arrival time of the signal = distance from
%         satellite to reciever + offset b
% OUTPUTS: eX estimated position of the user
%          eb estimate of the offset time b
%
% Implements the solution of Stephen Bancroft, "An Algebraic Solution of the GPS
% Equations," in IEEE Transactions on Aerospace and Electronic Systems,
% vol. AES-21, no. 1, pp. 56-59, Jan. 1985, doi: 10.1109/TAES.1985.310538.
% https://ieeexplore.ieee.org/document/4104017
%
% pdf at
% https://cas.tudelft.nl/Education/courses/ee4c03/assignments/indoor_loc/Bancroft85loc.pdf
%  Aaron Becker, Feb 23, 2021

%% Initialization
if nargin < 1
    % Number of satellites. The node positions are randomized. The transmitter
    % is randomly positioned outside of the area covered by the nodes.
    nSat = 4; % number of nodes
    dim = 3;  %TODO: failes with dim = 1.  Not sure why.  Wors for 2 and 3
    satPos = 50*rand(nSat,dim)-25;  % position of nodes
    if dim == 2
        th = 2*pi*rand; %random angle
        Xpos = (50+(50*rand)).*[cos(th) sin(th)];%user location (unknown to satellites)
    else
        Xpos = 150*rand(1,dim)-75;
    end
    b = 15+15*rand;  %  timing offset  (unknown to satellites)
    noiseLvl = 0; % Ideal signal + Noise.
    tVals = sqrt(sum((Xpos-satPos).^2,2))+b + noiseLvl*rand(nSat,1);  % Eqn (2)
else
    nSat = length(tVals);
    b = NaN;  % unknown to satellites
    Xpos = NaN(1,size(satPos,2)); % unknown to satellites
end

if nSat < size(satPos,2) +1
    disp(['Error: system is underdetermined with ',num2str(nSat),' satellites and ', num2str(size(satPos,2)), ' dimensions'])
end
%% calculation of user location and time offset
A = [satPos, tVals];                         % Eqn (5)
i0 = ones(nSat,1);                           % Eqn (6)
r = ones(nSat,1);                            % Eqn (7)
for i = 1:nSat
    r(i) = minkowskiSum(A(i,:),A(i,:))/2;   % Eqn (8)
end
B = pinv(A);  % same as inv(A'*A)*A';       % Eqn (9)
u = B*i0;                                   % Eqn (10)
v = B*r;                                    % Eqn (11)
E = minkowskiSum(u,u);                      % Eqn (12)
F = minkowskiSum(u,v) - 1;                  % Eqn (13)
G = minkowskiSum(v,v);                      % Eqn (14)

lam1 = (-2*F+sqrt((2*F)^2-4*E*G))/(2*E);    % Eqn (15) +solution to quadratic
lam2 = (-2*F-sqrt((2*F)^2-4*E*G))/(2*E);    % Eqn (15) -solution to quadratic
y1 = lam1*u+v;                              % Eqn (16) +
y2 = lam2*u+v;                              % Eqn (16) -

Tx1 = y1(1:end-1)';                         % Eqn (17)
Tx2 = y2(1:end-1)';                         % Eqn (17)

b1 = -y1(end);                              % Eqn (17)
b2 = -y2(end);                              % Eqn (17)

%determine which is right by calculating the average error
tvals1 = sqrt(sum((Tx1-satPos).^2,2))+b1;
tvals2 = sqrt(sum((Tx2-satPos).^2,2))+b2;

% e1 = mean(abs(tvals1-tVals));
% e2 = mean(abs(tvals2-tVals));
e1 = std(tvals1-tVals);
e2 = std(tvals2-tVals);

if e1 < e2
    legent = {'Satellite positions','X position','estimate 1 (correct)','estimate 2'};
    eX = Tx1;
    eb = b1;
else
    legent = {'Satellite positions','X position','estimate 1','estimate 2 (correct)'};
    eX = Tx2;
    eb = b2;
end

% %% Plot the process
% cplot = @(r,x0,y0,linet) plot(x0 + r*cos(linspace(0,2*pi)),y0 + r*sin(linspace(0,2*pi)),linet);
%  
% if size(satPos,2) == 1
%     figure(1); clf; hold on;
%     p(1) = plot(satPos,zeros(size(satPos)),'bx');
%     hold on
%     for i = 1:length(satPos)
%         cplot(tVals(i)-b,satPos(i),0,'g-');
%     end
%     p(2) = plot(Xpos(1),0, 'xg');
%     set(p(2),'color', [0,0.5,0])
%     xlabel('x');ylabel('y')
%     axis equal
%     if  abs(max(tVals)-min(tVals) - (max(satPos)-min(satPos))) < 1e-10
%         [minV,minInd] = min(satPos);
%         [maxV,maxInd] = max(satPos);
%         legent = {'Satellite positions','X position'};
%         if tVals(maxInd) < tVals(minInd)
%             str = sprintf('The receiver position is unknown, but > %.2f',maxV);
%             a = axis;
%             reg = patch([max(satPos),max(satPos),a(2),a(2)],[a(4),a(3),a(3),a(4)],[0.5,0.5,1]);
%             uistack(reg,'bottom')
%         else
%             str = sprintf('The receiver position is unknown, but < %.2f',minV);
%             a = axis;
%             reg = patch([min(satPos),min(satPos),a(1),a(1)],[a(4),a(3),a(3),a(4)],[0.5,0.5,1]);
%             uistack(reg,'bottom')
%         end
%         title(sprintf('Act X = [%.2f], b=%.2f \n %s',Xpos(1),b,str))
%     else
%         p(3) = plot(Tx1(1),0, 'ob');
%         p(4) = plot(Tx2(1),0, 'or'); 
%         for i = 1:length(satPos)
%             cplot(tVals(i)-b,satPos(i),0,'g-');
%             cplot(tVals(i)-b1,satPos(i),0,'b--');
%             cplot(tVals(i)-b2,satPos(i),0,'r--');
%         end
%         title(sprintf('Act X = [%.2f], b=%.2f \n Est X = [%.2f], b=%.2f',Xpos(1),b,eX(1),eb))
%     end
%     legend(p,legent);
%     axis tight
%     disp(['Act X = [',num2str(Xpos),'], b=',num2str(b)])
%     disp(['Est X = [',num2str(eX),  '], b=',num2str(eb)])
%     disp('std(error),   b,       Tx  for two solutions:')
%     disp([e1,b1,Tx1])
%     disp([e2,b2,Tx2])
% elseif size(satPos,2) == 2
%     figure(1); clf; hold on;
%     p(1) = plot(satPos(:,1),satPos(:,2),'bx');
%     hold on
%     for i = 1:length(satPos)
%         cplot(tVals(i)-b,satPos(i,1),satPos(i,2),'g-');
%         cplot(tVals(i)-b1,satPos(i,1),satPos(i,2),'b--');
%         cplot(tVals(i)-b2,satPos(i,1),satPos(i,2),'r--');
%     end
%     
%     p(2) = plot(Xpos(1),Xpos(2), 'x');
%     set(p(2),'color', [0,0.5,0])
%     p(3) = plot(Tx1(1),Tx1(2), 'ob');
%     p(4) = plot(Tx2(1),Tx2(2), 'or');
%     
%     legend(p,legent);
%     xlabel('x');ylabel('y')
%   
%     title(sprintf('Act X = [%.2f,%.2f], b=%.2f \n Est X = [%.2f,%.2f], b=%.2f',Xpos(1),Xpos(2),b,eX(1),eX(2),eb))
%     axis equal
% elseif size(satPos,2) == 3
%     figure(1); clf; hold on;
%     p(1) = plot3(satPos(:,1),satPos(:,2),satPos(:,3),'bx');
%     hold on
%     [X,Y,Z] = sphere;
%     lightGrey = 0.8*[1 1 1];
%     for i = 1:length(satPos)
%         radius = tVals(i)-b;
%         surf(satPos(i,1)+radius*X,satPos(i,2)+radius*Y,satPos(i,3)+radius*Z,'facecolor','g','facealpha',0.02,'EdgeColor',lightGrey,'edgealpha',0.2)
%     end
%     
%     p(2) = plot3(Xpos(1),Xpos(2),Xpos(3), 'x');
%     set(p(2),'color', [0,0.5,0])
%     p(3) = plot3(Tx1(1),Tx1(2),Tx1(3), 'ob');
%     p(4) = plot3(Tx2(1),Tx2(2),Tx2(3), 'or');
%     
%     legend(p,legent);
%     xlabel('x');ylabel('y');zlabel('z')
%   
%     title(sprintf('Act X = [%.2f,%.2f,%.2f], b=%.2f \n Est X = [%.2f,%.2f,%.2f], b=%.2f',Xpos(1),Xpos(2),Xpos(3),b,eX(1),eX(2),eX(3),eb))
%     axis equal
%     view([1 1 0.75]) % adjust the viewing angle
% elseif size(satPos,2) > 2
%     disp(['plot is only supported for 1D & 2D, this is ',num2str(size(satPos,2)),'D'])
%     disp(['Act X = [',num2str(Xpos),'], b=',num2str(b)])
%     disp(['Est X = [',num2str(eX),  '], b=',num2str(eb)])
%     disp('std(error),   b,       Tx  for two solutions:')
%     disp([e1,b1,Tx1])
%     disp([e2,b2,Tx2])
% end
    function mSum = minkowskiSum(a,b)
        mSum = sum(a(1:end-1).*b(1:end-1)) - a(end)*b(end);  % Eqn (4)
    end

end
