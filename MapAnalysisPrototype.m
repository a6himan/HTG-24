% Housekeeping
clear;
clc;

set(groot, 'DefaultAxesFontSize',14)

% Center point
latCenter = 2.9;
lonCenter = 31.5;

% Calculate corner points of the square
lat = [latCenter+0.0045, latCenter+0.0045, latCenter-0.0045, latCenter-0.0045, latCenter+0.0045];
lon = [lonCenter-0.0045, lonCenter+0.0045, lonCenter+0.0045, lonCenter-0.0045, lonCenter-0.0045];

% Initialize matrices to hold the latitudes and longitudes of the points
latPoints = linspace(latCenter-0.0045, latCenter+0.0045, 50);
lonPoints = linspace(lonCenter-0.0045, lonCenter+0.0045, 50);

[lonGrid, latGrid] = meshgrid(lonPoints, latPoints);
xCoordinates = reshape(lonGrid, [], 1)';
yCoordinates = reshape(latGrid, [], 1)';

[lonGrid, latGrid] = meshgrid(lonPoints, latPoints);
latPointsFlat = reshape(latGrid, [], 1);
lonPointsFlat = reshape(lonGrid, [], 1);


%Plot sattelite image

%Plot sattelite image
figure(1);
basemap = 'satellite';
geobasemap(basemap);

% Center point
%latCenter = 2.9;
%lonCenter = 31.5;

% Calculate corner points of the square
% lat = [latCenter+0.009, latCenter+0.009, latCenter-0.009, latCenter-0.009, latCenter+0.009];
% lon = [lonCenter-0.009, lonCenter+0.009, lonCenter+0.009, lonCenter-0.009, lonCenter-0.009];
% 
% Initialize matrices to hold the latitudes and longitudes of the points
% latPoints = linspace(latCenter-0.009, latCenter+0.009, 50);
% lonPoints = linspace(lonCenter-0.009, lonCenter+0.009, 50);
% 

geoplot(lat, lon, '-r', 'LineWidth', 2) % Plot the square with red lines
hold on;
geoscatter(latPointsFlat, lonPointsFlat, 20,'red','filled');
% Set geographic limits to focus on the square with some margin
geolimits([latCenter-0.005 latCenter+0.005], [lonCenter-0.005 lonCenter+0.005])
hold off;

% Add a title
title('1km x 1km Discretization region','FontSize',14)

gx = gca;
[latlim,lonlim] = geolimits(gx);

% Conduct elevation analysis
layer = wmsfind("mathworks",SearchField="serverurl");
layer = refine(layer,"elevation");
[Z,RZ] = wmsread(layer,Latlim=latlim,Lonlim=lonlim,ImageFormat="image/bil");
Z = double(Z);
heights = zeros(size(latPointsFlat));

% Calculate height at each point
for i = 1:length(latPointsFlat)
    heights(i) = geointerp(Z, RZ, latPointsFlat(i), lonPointsFlat(i), 'nearest');
end
% Reshape heights to a 10x10 matrix to match the grid
heightsMatrix = reshape(heights, [50, 50]);
heightsFlat = reshape(heightsMatrix,[],1);

% Concatenate the XYZ coordinates into a single matrix
GridCoordinates = [lonPointsFlat, latPointsFlat, heightsFlat];

% % Constants
% radius_m = 5; % Assuming you want a 3m radius for the circles
% EARTH_RADIUS_M = 6371000; % Average Earth radius in meters
% numPoints = 100;
% 
% 
% 
% figure(2);
% 
% % Xnorm = linspace(0,500,50);
% % Ynorm = linspace(0,500,50);
% % 
% % [X2,Y2] = meshgrid(Xnorm,Ynorm);
% % surf(X2,Y2,heightsMatrix);
surf(latGrid, lonGrid, heightsMatrix); % Use the grid matrices for X and Y
% hold on
% % 
% % for i = 1 : 90 : length(latPointsFlat)
% % % Convert radius from meters to degrees of latitude
% % radius_deg_lat = (radius_m / 63710000) * 180 / pi;
% % radius_deg_lon = radius_deg_lat / cosd(latPointsFlat(i));
% % 
% % % Generate the circle in degrees
% % theta = linspace(0, 2*pi, numPoints);
% % latCircle = latPointsFlat(i) + (radius_deg_lat * cos(theta)) * (180 / pi);
% % lonCircle = lonPointsFlat(i) + (radius_deg_lon * sin(theta)) * (180 / pi) / cosd(latPointsFlat(i));
% % 
% % % Interpolate z values for the circle perimeter directly
% % zCircle = interp2(lonGrid, latGrid, heightsMatrix, lonCircle, latCircle, 'cubic');
% % 
% % 
% % % Assuming the surf plot is the current figure
% % plot3(latCircle, lonCircle, zCircle, 'r-', 'LineWidth', 2);
% % end
% 
% hold off
xlabel('Latitude (decimal °)','FontSize',14);
ylabel('Longitude (decimal °)','FontSize',14);
zlabel('Elevation above sea level (m)','FontSize',14);
title('3D Elevation Visualisation','FontSize',14);
colorbar;

% 
% % Constants
% radius_m = 15; % Assuming you want a 3m radius for the circles
% EARTH_RADIUS_M = 6371000; % Average Earth radius in meters
% numPoints = 100;
% CircleLonCenter = lonCenter;
% CircleLatCenter = latCenter;
% 
% 
% % Convert radius from meters to degrees of latitude
% radius_deg_lat = (radius_m / 67310000) * 180 / pi;
% radius_deg_lon = radius_deg_lat / cosd(CircleLatCenter);
% 
% % Generate the circle in degrees
% theta = linspace(0, 2*pi, numPoints);
% latCircle = CircleLatCenter + (radius_deg_lat * cos(theta)) * (180 / pi);
% lonCircle = CircleLonCenter + (radius_deg_lon * sin(theta)) * (180 / pi) / cosd(CircleLatCenter);
% 
% % Interpolate z values for the circle perimeter directly
% zCircle = interp2(lonGrid, latGrid, heightsMatrix, lonCircle, latCircle, 'cubic');
% 
% figure(2);
% %surf(latGrid, lonGrid, heightsMatrix); % Use the grid matrices for X and Y
% %hold on
% 
% % Assuming the surf plot is the current figure
% plot3(latCircle, lonCircle, zCircle, 'r-', 'LineWidth', 2);
% %hold off
% xlabel('Latitude (decimal °)');
% ylabel('Longitude (decimal °)');
% zlabel('Elevation above sea level (m)');
% title('3D water sprinkler range');
% grid on; box on; 
% %colorbar;
% %figure(3)
% %plot3(latCircle, lonCircle, zCircle, 'r-', 'LineWidth', 2);





