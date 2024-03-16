function [latCircle, lonCircle, zCircle] = generateCirclePerimeter(CircleLatCenter, CircleLonCenter, GridCoordinates)
    
% Constants
radius_m = 30; % Assuming you want a 3m radius for the circles
EARTH_RADIUS_M = 6371000; % Average Earth radius in meters
numPoints = 100;
%CircleLonCenter = lonCenter;
%CircleLatCenter = latCenter;


% Convert radius from meters to degrees of latitude
radius_deg_lat = (radius_m / 67310000) * 180 / pi;
radius_deg_lon = radius_deg_lat / cosd(CircleLatCenter);

% Generate the circle in degrees
theta = linspace(0, 2*pi, numPoints);
latCircle = CircleLatCenter + (radius_deg_lat * cos(theta)) * (180 / pi);
lonCircle = CircleLonCenter + (radius_deg_lon * sin(theta)) * (180 / pi) / cosd(CircleLatCenter);

% Interpolate z values for the circle perimeter directly
zCircle = interp2(lonGrid, latGrid, heightsMatrix, lonCircle, latCircle, 'cubic');
end

