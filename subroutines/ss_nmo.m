function [dout,t410_ref,t660_ref] = ss_nmo(din,t,h,h0,is_plot)
% move-out correction for SS precursors
% Input:
% din: nt*1 vector of a SS trace
% t: nt*1 vector of time axis
% h: epicenter distance of SS trace
% h0: reference distance needs to be corrected to
% is_plot: flag that controls plotting
% Output:
% dout: nt*1 vector that contains NMO corrected trace
% Jan. 21, 2022, Yunfeng Chen, Global Seismology Group, Zhejiang University

% phase arrival time of targeting trace
times=taupTime('ak135',10,'SS,S^410S,S^660S','deg',h);
indices = find(strcmp({times.phaseName},'S^660S'));
t660 =times(indices(1)).time;
indices = find(strcmp({times.phaseName},'S^410S'));
t410 =times(indices(1)).time;
indices = find(strcmp({times.phaseName},'SS'));
tss =times(indices(1)).time;
t410=t410-tss;
t660=t660-tss;
t520=(t660+t410)/2;

% phase arrival time of reference trace
times=taupTime('ak135',10,'SS,S^410S,S^660S','deg',h0);
indices = find(strcmp({times.phaseName},'S^660S'));
t660_ref =times(indices(1)).time;
indices = find(strcmp({times.phaseName},'S^410S'));
t410_ref =times(indices(1)).time;
indices = find(strcmp({times.phaseName},'SS'));
tss_ref =times(indices(1)).time;
t410_ref=t410_ref-tss_ref;
t660_ref=t660_ref-tss_ref;
dt=t(2)-t(1);

npoints=30;
dseg1=din(t<t520+npoints);
tseg1=t(t<t520+npoints);
dseg2=din(t>=t520-npoints);
tseg2=t(t>=t520-npoints);
tc1=t660_ref-t660;
dseg1_shift = fftShift(dseg1,tseg1,tc1);
tc2=t410_ref-t410;
dseg2_shift = fftShift(dseg2,tseg2,tc2);
d_shift=[dseg1_shift(1:end-npoints); dseg2_shift(npoints+1:end)];
dout=d_shift;
if is_plot
    figure;
    plot(t,din,t,dout);
end
