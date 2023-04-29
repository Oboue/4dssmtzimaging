% Jun. 11, 2022, Yunfeng Chen, Global Seismology Group, Zhejiang University
% crustal correction of ss precursor
function [corr,ttmod,ttref] = ss_crucor(lat,lon,rayp)
elev=topo_points(lat,lon,'crust1.0');
elev=elev/1000;
[corr,ttmod,ttref]=crucor(lat,lon,rayp,'S','refmod','prem','elev',elev);
corr=corr*2;
ttmod=ttmod*2;
ttref=ttref*2;
% plot(tcorr,dcdps,'.'); hold on;
% plot(tc,z)