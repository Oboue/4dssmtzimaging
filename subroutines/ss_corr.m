% Jun. 6, 2022, Yunfeng Chen, Global Seismology Group, Zhejiang University
% heterogeneity correction of ss precursor
function [tc,vs3d] = ss_corr(h0, lon0, lat0, Fvs, z, vs, zi)
zmid=diff(zi)/2+zi(1:end-1);
vs1d = interp1(z,vs,zmid);
% calcualte the raypath
raypath=taupPath('prem',0,'SS,S^400S','deg',h0);
xp=raypath.path.distance;
zp=raypath.path.depth;
ztmp=zp(zp<=1000 & xp>=3*h0/4);
xtmp=xp(zp<=1000 & xp>=3*h0/4);
xi=interp1(ztmp,xtmp,zi);
dl=sqrt(diff(zi).^2+diff(xi).^2);
% extract 1D velocity model
xi=lon0-2.5:0.5:lon0+2.5;
yi=lat0-2.5:0.5:lat0+2.5;
[XI,YI,ZI]=meshgrid(xi,yi,zmid);
VI=Fvs(XI,YI,ZI);
vs3d=squeeze(mean(mean(VI,1),2));
% figure;
% for ii=1:5
%     for jj=1:5
%         plot(squeeze(VI(ii,jj,:))); hold on;
%     end
% end
% plot(dvs,'linewidth',3)

idx=isnan(vs3d);
vs3d(idx)=vs1d(idx);

%         figure;
%         plot(zmid,vs1d,zmid,vs3d);

tss3d = dl./vs3d;
tss1d = dl./vs1d;

temp = tss1d-tss3d;
temp(isnan(temp)) = 0;
tc = cumsum([0; temp]);
%         figure;
%         plot(zi,TimeCorrections);