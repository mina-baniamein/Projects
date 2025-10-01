function plotearth()
[X,Y,Z]=sphere(50);
R=6378;
globe= surf(-X*R,Y*R,-Z*R);
cdata = imread('earth.jpg');
set(globe, 'FaceColor', 'texturemap', 'CData', cdata,  'EdgeColor', 'none');
set(gcf,'Color','w')
set(gca, 'visible', 'on')
axis equal
view(-29,9)
end