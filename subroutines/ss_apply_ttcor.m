function [dout] = ss_apply_ttcor(din,t,tt1d,tc,is_plot)
% apply travel time correction
% Input:
% din: nt*1 vector of a SS trace
% t: nt*1 vector of time axis
% tt1d: differential travel time in 1D velocity model TSS1D-TSdS1D
% tt3d: differential travel time in 3D velocity model TSS3D-TSdS3D
% is_plot: flag that controls plotting
% Output:
% dout: nt*1 vector that contains NMO corrected trace
% Jun. 11, 2022, Yunfeng Chen, Global Seismology Group, Zhejiang University
% phase arrival time of targeting trace
if nargin<5
    is_plot=false;
end
t=t(:);
dtmp=interp1(t,din,tt1d);
dtmp=[dtmp(end:-1:1), din(t>0)'];
% interpolate again to obtain regular time axis
tt3d=tt1d+tc;
tnew=t;
t0=[tt3d(end:-1:1), t(t>0)']';
dout=interp1(t0,dtmp,tnew,'linear',0);
if is_plot
    figure;
    plot(t,din,t,dout);
end