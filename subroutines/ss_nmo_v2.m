function [dout,t410_ref,t660_ref] = ss_nmo_v2(din,t,h,h0,F,is_plot)
% move-out correction for SS precursors using travel time table to stretch
% and squeeze the trace
% Input:
% din: nt*1 vector of a SS trace
% t: nt*1 vector of time axis
% h: epicenter distance of SS trace
% h0: reference distance needs to be corrected to
% F: interpolation function (nz*nx) that contains the travel time for each
%    depth and distance
% is_plot: flag that controls plotting
% Output:
% dout: nt*1 vector that contains NMO corrected trace
% Jan. 23, 2022, Yunfeng Chen, Global Seismology Group, Zhejiang University
% phase arrival time of reference trace
y0 = 0:1:1000; 
x0=h0*ones(size(y0));
t=t(:);
t0=F(x0,y0);
t0=[t0(end:-1:1),t(t>0)'];
% phase arrival time of targeting trace
y=y0;
x=h*ones(size(y));
ti=F(x,y);
dtmp=interp1(t,din,ti);
dtmp=[dtmp(end:-1:1), din(t>0)'];
% interpolate again to obtain regular time axis
tnew=t;
dout=interp1(t0,dtmp,tnew,'linear',0);
t410_ref=F(h0,410);
t660_ref=F(h0,660);
if is_plot
    figure;
    plot(t,din,t,dout);
end
