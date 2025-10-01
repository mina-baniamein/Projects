%% Textured 3D Earth example
%
% Ryan Gray
% 8 Sep 2004
% Revised 9 March 2006, 31 Jan 2006, 16 Oct 2013
%% Options
space_color = 'w';
npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible
%GMST0 = []; % Don't set up rotatable globe (ECEF)
GMST0 = 4.89496121282306; % Set up a rotatable globe at J2000.0
% Earth texture image
% Anything imread() will handle, but needs to be a 2:1 unprojected globe
% image.
image_file = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';
% Mean spherical earth
erad    = 6371.0087714; % equatorial radius (meters)
prad    = 6371.0087714; % polar radius (meters)
erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)
%% Create figure
figure('Color', space_color);
hold on;
% Turn off the normal axes
%set(gca, 'NextPlot','add', 'Visible','off');
%axis equal;
%axis auto;
% Set initial view
view(0,30);
axis vis3d;
%% Create wireframe globe
% Create a 3D meshgrid of the sphere points using the ellipsoid function
[x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
hold on
% for j=1:lenght(rr)
%  plot3(rrplot(1,:),rrplot(2,:),rrplot(3,:),'b','lineWidth',2);
%  plot3(rrplot_f(1,:),rrplot_f(2,:),rrplot_f(3,:),'r','lineWidth',2)
%  legend('orbita iniziale','orbita finale')
% pause(0.1);
% %comet3(vvplot(1,:),vvplot(2,:),vvplot(3,:));
% end
globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
if ~isempty(GMST0)
    hgx = hgtransform;
    set(hgx,'Matrix', makehgtform('zrotate',GMST0));
    set(globe,'Parent',hgx);
end
%% Texturemap the globe
% Load Earth image for texture map
cdata = imread(image_file);
% Set image as color data (cdata) property, and set face color to indicate
% a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');

% for j=1:length(rrplot)
% plot3(rrplot(1,j),rrplot(2,j),rrplot(3,j));
% hold on
% pause(0.001);
% %comet3(vvplot(1,:),vvplot(2,:),vvplot(3,:));
% end
