function [bplat,bplon]=gc_midpoint(lat1, lon1, lat2, lon2)
lat1r = deg2rad(lat1);
lat2r = deg2rad(lat2);
lon1r = deg2rad(lon1);
lon2r = deg2rad(lon2);

bx = cos(lat2r) * cos(lon2r-lon1r);
by = cos(lat2r) * sin(lon2r-lon1r);

bplatr = atan2( sin(lat1r)+sin(lat2r), sqrt((cos(lat1r)+bx)^2+by^2 ));
bplonr = lon1r + atan2(by, cos(lat1r)+bx);

bplat = rad2deg(bplatr);
bplon = rad2deg(bplonr);

