%{
  Use this to test convex_hull_3D
%}

%%%%%%%%%%%%%%%%%%
%%%% 3D FANS  %%%%
%%%%%%%%%%%%%%%%%%

% CHOOSE SET OF POINTS:
%{
% CUBIC (ALWAYS A CONVEX SET)
p3 = [0 0 0;
      0 0 1;
      0 1 0;
      0 1 1;
      1 0 0;
      1 0 1;
      1 1 0;
      1 1 1];
% WITH RANDOM VARIANCE
%p3 = p3 + 0.5*rand(size(p3));
%}
% SPHERICAL (ALWAYS A CONVEX SET)
n = 20;
a = [2*pi*rand(n,1) pi*rand(n,1)];
p3 = [cos(a(:,1)).*sin(a(:,2)) sin(a(:,1)).*sin(a(:,2)) cos(a(:,2))];
%{
% RANDOM POINTS (NEVER A CONVEX SET)
n = 20;
p3 = rand(n,3);
% TO GUARANTEE NOT CONVEX SET
p3(n,:) = mean(p3(1:n-1,:));
%}
%

% CALCULATIONS
l3 = convex_hull_3D(p3);

%
figure(1)
clf

% CHOOSE PLOT STYLE
%{
% CONVEX HULL (VANILLA)
colour = 0.5*ones(size(p3));
trisurf(l3, p3(:,1), p3(:,2), p3(:,3), colour)
%}
% SEPARATED PIECES (TO SEE THAT THERE ARE NO OVERLAPS OR PIECES INSIDE)
hold on
prec = 0.05;
for ii = 1:size(p3,1)
	l3s = l3(find(l3(:,1)==ii),:);
	% colour
	colour = ii/size(p3,1)*ones(size(p3));
	% plot
	trisurf(l3s, p3(:,1)+prec*p3(ii,1), p3(:,2)+prec*p3(ii,2), p3(:,3)+prec*p3(ii,3), colour)
end
hold off
%

box off
%

% AUTOROTATION
%
secs = 2; % period of rotation, edit this
rounds = 1; % number of revolutuions
msecs = 0.05;
degs = msecs/secs * 360;
eldegs = 45;
az = 0;
el = 0;
view(az,el)
while az < 360*rounds
	pause(msecs)
	az = az + degs;
	el = eldegs * sin(az/41);
	view(az,el)
end
%
