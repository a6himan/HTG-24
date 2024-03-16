% Housekeeping
clear;
clc;

% Center point
latCenter = 2.9;
lonCenter = 31.5;

% Calculate corner points of the square
lat = [latCenter+0.009, latCenter+0.009, latCenter-0.009, latCenter-0.009, latCenter+0.009];
lon = [lonCenter-0.009, lonCenter+0.009, lonCenter+0.009, lonCenter-0.009, lonCenter-0.009];

% Initialize matrices to hold the latitudes and longitudes of the points
latPoints = linspace(latCenter-0.009, latCenter+0.009, 10);
lonPoints = linspace(lonCenter-0.009, lonCenter+0.009, 10);

[lonGrid, latGrid] = meshgrid(lonPoints, latPoints);
xCoordinates = reshape(lonGrid, [], 1)';
yCoordinates = reshape(latGrid, [], 1)';


%Plot sattelite image

%Plot sattelite image
figure(1);
basemap = 'satellite';
geobasemap(basemap);

% Center point
%latCenter = 2.9;
%lonCenter = 31.5;

% Calculate corner points of the square
lat = [latCenter+0.009, latCenter+0.009, latCenter-0.009, latCenter-0.009, latCenter+0.009];
lon = [lonCenter-0.009, lonCenter+0.009, lonCenter+0.009, lonCenter-0.009, lonCenter-0.009];

% Initialize matrices to hold the latitudes and longitudes of the points
latPoints = linspace(latCenter-0.009, latCenter+0.009, 50);
lonPoints = linspace(lonCenter-0.009, lonCenter+0.009, 50);

[lonGrid, latGrid] = meshgrid(lonPoints, latPoints);
latPointsFlat = reshape(latGrid, [], 1);
lonPointsFlat = reshape(lonGrid, [], 1);
geoplot(lat, lon, '-r', 'LineWidth', 2) % Plot the square with red lines
hold on;
geoscatter(latPointsFlat, lonPointsFlat, 20,'red','filled');
% Set geographic limits to focus on the square with some margin
geolimits([latCenter-0.01 latCenter+0.01], [lonCenter-0.01 lonCenter+0.01])

% Add a title
title('1 km by 1 km Square in Northwest Uganda')

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



% Constants
radius_m = 30; % Assuming you want a 3m radius for the circles
EARTH_RADIUS_M = 6371000; % Average Earth radius in meters
numPoints = 100;
CircleLonCenter = lonCenter;
CircleLatCenter = latCenter;


% Convert radius from meters to degrees of latitude
radius_deg_lat = (radius_m / 67310000) * 180 / pi;
radius_deg_lon = radius_deg_lat / cosd(CircleLatCenter);

% Generate the circle in degrees
theta = linspace(0, 2*pi, numPoints);
latCircle = CircleLatCenter + (radius_deg_lat * cos(theta)) * (180 / pi);
lonCircle = CircleLonCenter + (radius_deg_lon * sin(theta)) * (180 / pi) / cosd(CircleLatCenter);

% Interpolate z values for the circle perimeter directly
zCircle = interp2(lonGrid, latGrid, heightsMatrix, lonCircle, latCircle, 'cubic');

figure(2);
surf(latGrid, lonGrid, heightsMatrix); % Use the grid matrices for X and Y
hold on

% Assuming the surf plot is the current figure
plot3(latCircle, lonCircle, zCircle, 'r-', 'LineWidth', 2);
hold off
xlabel('Latitude');
ylabel('Longtude');
zlabel('Elevation');
title('Elevation Surface Plot');
colorbar;

%figure;
%plot3(latCircle, lonCircle, zCircle, 'r-', 'LineWidth', 2);




