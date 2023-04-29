% Feb. 16, 2023, Yunfeng Chen, Global Seismology Group, Zhejiang University
% read in synthetic data from Maria Koroni and Jeannot Trampert.
% 
% Reference: Oboue Y.A.S.I, Y. Chen, J. Wang, X. Jiang, Ramin M.H. Dokht, Y. J. Gu, M. Koroni, and
% Y. Chen, 2023, High-resolution mantle transition zone imaging using
% multi-dimensional reconstruction of SS precursors, JGR, submitted

clear; close all; clc;
%% Read in data and calculate bounce point (midpoint)
addpath('../rdrr')
addpath('../data') 
javaaddpath('./FMI/lib/FMI.jar');
addpath('~/MATLAB/m_map');
addpath('../slab')
addpath('../MatSAC');
addpath('../utils');
addpath('../subroutines');
addpath('./FMI/matTaup');
addpath('../TX2019slab')
addpath export_fig

%% 4D pre-stack reconstruction using RDRR algorithm
disp('4D RDRR')

% load 4D synthetic data

load synthssp4d.mat 

% add random noise
randn('state',201314);

var=0.1;
synthssp4dn=d+var*randn(size(d));

% Remove 70 percent of traces 
ratio=0.3; % 70 percent missing traces 
mask=genmask(reshape(synthssp4dn,nt,nx*ny*nh),ratio,'c',201415);
mask=reshape(mask,nt,nx,ny,nh);
% size(ss4din)
synthssp4d0=synthssp4dn.*mask;

dt=t(2)-t(1);
flow=1/75;
fhigh=1/15.;
Niter=10;
mode=1;
verb=1;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing

N=33; % Rank 
K=2;  % Damping factor default 2 larger value preserve more signal 
u=2.7; % Cooling factor level of sparsity can be fixed
e=0.954;  % Rational transfer function coefficient can be fixed
ws=1;    % Window size 1 by default
iflb=0;
dsynthrec4d=rdrr5d_lb_recon(synthssp4d0,mask,flow,fhigh,dt,N,K,Niter,eps,verb,mode,iflb,a,u,e,ws);
% save synt4drecons.mat dsynthrec4d

o=yc_snr(d,synthssp4d0,2)
o=yc_snr(d,dsynthrec4d,2)

figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w'); subplot(4,1,1);
imagesc(1:nx*ny*nh,t,reshape(d,nt,nx*ny*nh)); hold on;
plot(1:nx*ny*nh,t410_ref*ones(1,nx*ny*nh),'--r','linewidth',4)
plot(1:nx*ny*nh,t660_ref*ones(1,nx*ny*nh),'--r','linewidth',4)
title('Pre-stack (true)')
colormap(seis);caxis([-0.12 0.12]); colorbar;
axis xy
ylim([-300 -100])
ylabel('Time to SS (sec)')
xlabel('Trace');
set(gca,'fontsize',16)
annotation(gcf,'textbox',...
    [0.781666666666666 0.880677412772908 0.073888888888889 0.0246548323471352],...
    'String','S410S',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
annotation(gcf,'textbox',...
    [0.783333333333334 0.828402366863905 0.0722222222222222 0.0246548323471401],...
    'String','S660S',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
text(-325,-75,'(a)','color','k','Fontsize',22,'fontweight','normal','HorizontalAlignment','center');

% figure;
subplot(4,1,2);
imagesc(1:nx*ny*nh,t,reshape(synthssp4d0,nt,nx*ny*nh)); hold on;
plot(1:nx*ny*nh,t410_ref*ones(1,nx*ny*nh),'--r','linewidth',3)
plot(1:nx*ny*nh,t660_ref*ones(1,nx*ny*nh),'--r','linewidth',3)
title('Pre-stack (raw)')
colormap(seis);caxis([-0.12 0.12]); colorbar;
axis xy
ylim([-300 -100])
ylabel('Time to SS (sec)')
xlabel('Trace');
set(gca,'fontsize',16)
annotation(gcf,'textbox',...
    [0.780555555555555 0.66272869482419 0.073888888888889 0.0246548323471352],...
    'String','S410S',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
annotation(gcf,'textbox',...
    [0.783333333333334 0.606511790596674 0.0722222222222222 0.0246548323471401],...
    'String','S660S',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
text(-325,-75,'(b)','color','k','Fontsize',22,'fontweight','normal','HorizontalAlignment','center');

subplot(4,1,3);
imagesc(1:nx*ny*nh,t,reshape(dsynthrec4d,nt,nx*ny*nh)); hold on;
plot(1:nx*ny*nh,t410_ref*ones(1,nx*ny*nh),'--r','linewidth',4)
plot(1:nx*ny*nh,t660_ref*ones(1,nx*ny*nh),'--r','linewidth',4)
title('Pre-stack (reconstructed)')
colormap(seis);caxis([-0.12 0.12]); colorbar;
axis xy
ylim([-300 -100])
ylabel('Time to SS (sec)')
xlabel('Trace');
set(gca,'fontsize',16)
annotation(gcf,'textbox',...
    [0.78 0.442802732192656 0.0755555555555556 0.0256410256410315],...
    'String','S410S',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
annotation(gcf,'textbox',...
    [0.783333333333333 0.394478289173251 0.0722222222222222 0.02465483234714],...
    'String','S660S',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
text(-325,-75,'(c)','color','k','Fontsize',22,'fontweight','normal','HorizontalAlignment','center');

error=d-dsynthrec4d;

subplot(4,1,4);
imagesc(1:nx*ny*nh,t,reshape(error,nt,nx*ny*nh)); hold on;
plot(1:nx*ny*nh,t410_ref*ones(1,nx*ny*nh),'--r','linewidth',4)
plot(1:nx*ny*nh,t660_ref*ones(1,nx*ny*nh),'--r','linewidth',4)
title('Reconstruction error')
colormap(seis);caxis([-0.12 0.12]); colorbar;
axis xy
ylim([-300 -100])
ylabel('Time to SS (sec)')
xlabel('Trace');
set(gca,'fontsize',16)
annotation(gcf,'textbox',...
    [0.78 0.22287968441815 0.0755555555555556 0.0256410256410316],...
    'String','S410S',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
annotation(gcf,'textbox',...
    [0.782222222222222 0.167652859960554 0.0722222222222222 0.02465483234714],...
    'String','S660S',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
text(-325,-75,'(d)','color','k','Fontsize',22,'fontweight','normal','HorizontalAlignment','center');

% o=yc_snr(d,synthssp4dn,2)

% [d4d_recon] = reconstruction(d4d,mask,dt,flow,fhigh,N,0.7,Niter);
% save output/d4dfield_recon_srn3_r35K3u0001e080_2.5deg.mat d4d d4d_recon x y h t 
%% stacking process
% post-stack (true)

d4d_stack = zeros(nt,nx,ny);
for i=1:nx
    for j=1:ny
        % move-out corrected cmp
        d_cmp_obs = squeeze(d(:,i,j,:));
        nstack = sum(any(d_cmp_obs));
        if nstack>0
            d_obs_stack = sum(d_cmp_obs,2)/nstack;
            d4d_stack(:,i,j)=d_obs_stack;
        end
    end 
end

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w'); 
subplot(4,1,1);
imagesc(1:nx*ny,t,reshape(d4d_stack,nt,nx*ny)); hold on;
plot(1:nx*ny,t410_ref*ones(1,nx*ny),'--r','linewidth',4)
plot(1:nx*ny,t660_ref*ones(1,nx*ny),'--r','linewidth',4)
title('Post-stack (true)')
colormap(seis);caxis([-0.12 0.12]); colorbar;
axis xy
ylim([-300 -100])
ylabel('Time to SS (sec)')
xlabel('Trace');
set(gca,'fontsize',16)
annotation(gcf,'textbox',...
    [0.777222222222222 0.883629191321505 0.073888888888889 0.0256410256410314],...
    'String','S410S',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
annotation(gcf,'textbox',...
    [0.777777777777778 0.828402366863906 0.0722222222222222 0.0246548323471401],...
    'String','S660S',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
text(-10,-75,'(a)','color','k','Fontsize',22,'fontweight','normal','HorizontalAlignment','center');

% save d4d_stack.mat d4d_stack t410_ref t660_ref nx ny nh t 

% post-stack (raw)
d4draw_stack = zeros(nt,nx,ny);
for i=1:nx
    for j=1:ny
        % move-out corrected cmp
        d_cmp_obs = squeeze(synthssp4d0(:,i,j,:));
        nstack = sum(any(d_cmp_obs));
        if nstack>0
            d_true_stack = sum(d_cmp_obs,2)/nstack;
            d4draw_stack(:,i,j)=d_true_stack;
        end
    end 
end

subplot(4,1,2);
% set(gcf,'Position',[100 100 1600 500],'color','w')
imagesc(1:nx*ny,t,reshape(d4draw_stack,nt,nx*ny)); hold on;
plot(1:nx*ny,t410_ref*ones(1,nx*ny),'--r','linewidth',3)
plot(1:nx*ny,t660_ref*ones(1,nx*ny),'--r','linewidth',3)
title('Post-stack (raw)')
colormap(seis);caxis([-0.12 0.12]); colorbar;
axis xy
ylim([-300 -100])
ylabel('Time to SS (sec)')
xlabel('Trace');
set(gca,'fontsize',16)
annotation(gcf,'textbox',...
    [0.777222222222222 0.669625246548328 0.073888888888889 0.0246548323471352],...
    'String','S410S',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
annotation(gcf,'textbox',...
    [0.777777777777778 0.609467455621305 0.0722222222222222 0.0246548323471401],...
    'String','S660S',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
text(-10,-75,'(b)','color','k','Fontsize',22,'fontweight','normal','HorizontalAlignment','center');

% save d4draw_stack.mat d4draw_stack t410_ref t660_ref nx ny nh t 

% post-stack (reconstructed)
d4d_recon_stack = zeros(nt,nx,ny);
for i=1:nx
    for j=1:ny
%         move-out corrected cmp
        d_cmp = squeeze(dsynthrec4d(:,i,j,:));
        nstack = sum(any(d_cmp));
        if nstack>0
            d4d_rec_stack= sum(d_cmp,2)/nstack;
            d4d_recon_stack(:,i,j)=d4d_rec_stack;
        end
    end 
end

subplot(4,1,3);
imagesc(1:nx*ny,t,reshape(d4d_recon_stack,nt,nx*ny)); hold on;
plot(1:nx*ny,t410_ref*ones(1,nx*ny),'--r','linewidth',4)
plot(1:nx*ny,t660_ref*ones(1,nx*ny),'--r','linewidth',4)
title('Post-stack (reconstructed)')
colormap(seis);caxis([-0.12 0.12]); colorbar;
axis xy
ylim([-300 -100])
ylabel('Time to SS (sec)')
xlabel('Trace');
set(gca,'fontsize',16)
annotation(gcf,'textbox',...
    [0.775555555555556 0.44477317554241 0.0755555555555556 0.0256410256410315],...
    'String','S410S',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
annotation(gcf,'textbox',...
    [0.778888888888889 0.388560157790929 0.0722222222222222 0.02465483234714],...
    'String','S660S',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
text(-10,-75,'(c)','color','k','Fontsize',22,'fontweight','normal','HorizontalAlignment','center');

% save d4d_recon_stack.mat d4d_recon_stack t410_ref t660_ref nx ny nh t 
% post-stack --- reconstruction error

error=d4d_stack-d4d_recon_stack;

subplot(4,1,4);
imagesc(1:nx*ny,t,reshape(error,nt,nx*ny)); hold on;
plot(1:nx*ny,t410_ref*ones(1,nx*ny),'--r','linewidth',4)
plot(1:nx*ny,t660_ref*ones(1,nx*ny),'--r','linewidth',4)
title('Reconstruction error')
colormap(seis);caxis([-0.12 0.12]); colorbar;
axis xy
ylim([-300 -100])
ylabel('Time to SS (sec)')
xlabel('Trace');
set(gca,'fontsize',16)
annotation(gcf,'textbox',...
    [0.774444444444444 0.22287968441815 0.0755555555555556 0.0256410256410316],...
    'String','S410S',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
annotation(gcf,'textbox',...
    [0.777777777777778 0.167652859960554 0.0722222222222222 0.02465483234714],...
    'String','S660S',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
text(-10,-75,'(d)','color','k','Fontsize',22,'fontweight','normal','HorizontalAlignment','center');

o=yc_snr(d4d_stack,d4draw_stack,2)
o=yc_snr(d4d_stack,d4d_recon_stack,2)
% save output/d4dfield_recon_stack_srn3_r35K3u0001e080_2.5deg.mat d4d_stack d4d_recon_stack x y t 
%% plot figures pre-stack raw and reconstructed 

% size(d4d_recon)
% time to depth conversion
dist = 95:5:170;
depth = 0:5:1000; 
[tt, f]=ss_tt_table(dist,depth);
z0=0:1:1000; 
h0=135;
x0=h0*ones(size(z0));
t0=f(x0,z0);

keep = t>=-300 & t<=-100;

t=t(keep);
ti=t;
% t0=-300;
% t1=-100;
% %t0=-500;
% %t1=100;
% keep = t>=t0 & t<=t1;

nz=length(z0);
d3d0_depth=zeros(nz,nx,ny);
d3d1_depth=zeros(nz,nx,ny);
d3d2_depth=zeros(nz,nx,ny);
% d3d3_depth=zeros(nz,nx,ny);

for i = 1:nx
    for j = 1:ny
        d=d4d_stack(:,i,j);
        dtmp=interp1(ti,d,t0,'linear',0);
        d3d0_depth(:,i,j)=dtmp;
        
        d=d4draw_stack(:,i,j);
        dtmp=interp1(ti,d,t0,'linear',0);
        d3d1_depth(:,i,j)=dtmp;
        
        d=d4d_recon_stack(:,i,j);
        dtmp=interp1(ti,d,t0,'linear',0);
        d3d2_depth(:,i,j)=dtmp;

%         d=d3(:,i,j);
%         dtmp=interp1(ti,d,t0,'linear',0);
%         d3d3_depth(:,i,j)=dtmp;
    end
end

%plot the trace in depth domain
d2d0_depth=reshape(d3d0_depth,nz,nx*ny);
d2d1_depth=reshape(d3d1_depth,nz,nx*ny);
d2d2_depth=reshape(d3d2_depth,nz,nx*ny);

a=yc_snr(d3d0_depth,d3d1_depth,2)

b=yc_snr(d3d0_depth,d3d2_depth,2)

% A = rms(d2d0_depth,d2d1_depth)
% B = rms(d2d0_depth,d2d2_depth)

% d2d3_depth=reshape(d3d3_depth,nz,nx*ny);
% figure;
% set(gcf,'Position',[0 0 1200 1000],'Color','w')
% subplot(411)
% imagesc(x,z0,d2d0_depth)
% ylim([200 800])
% % axis xy;
% caxis([-0.08 0.08])
% xlabel('CMP#')
% ylabel('Depth (km)')
% % title('r=150, snr<=3')
% set(gca,'fontsize',14)
% subplot(412)
% imagesc(x,z0,d2d1_depth)
% ylim([200 800])
% caxis([-0.08 0.08])
% xlabel('CMP#')
% ylabel('Depth (km)')
% % title('r=150, snr<=0.5')
% set(gca,'fontsize',14)
% subplot(413)
% imagesc(x,z0,d2d2_depth)
% ylim([200 800])
% colormap(seismic(3))
% caxis([-0.08 0.08])
% xlabel('CMP#')
% ylabel('Depth (km)')
% % title('r=500, snr<=3')
% set(gca,'fontsize',14)
% subplot(414)
% imagesc(x,z0,d2d3_depth)
% ylim([200 800])
% colormap(seismic(3))
% caxis([-0.08 0.08])
% xlabel('CMP#')
% ylabel('Depth (km)')
% % title('r=500, snr<=0.5')
% set(gca,'fontsize',14)
% export_fig 'post-stack_depth.png' '-r150'

%% find the maximum amplitude in the given depth intervals
depth0=z0;
drange1=[390,430];
drange2=[640,680];
drange3=[500,540];
% 
amp410_d0=zeros(nx,ny);
amp520_d0=zeros(nx,ny);
amp660_d0=zeros(nx,ny);
d410_d0=zeros(nx,ny);
d520_d0=zeros(nx,ny);
d660_d0=zeros(nx,ny);

amp410_d1=zeros(nx,ny);
amp520_d1=zeros(nx,ny);
amp660_d1=zeros(nx,ny);
d410_d1=zeros(nx,ny);
d520_d1=zeros(nx,ny);
d660_d1=zeros(nx,ny);
% 
amp410_d2=zeros(nx,ny);
amp520_d2=zeros(nx,ny);
amp660_d2=zeros(nx,ny);
d410_d2=zeros(nx,ny);
d520_d2=zeros(nx,ny);
d660_d2=zeros(nx,ny);
% 
amp410_d3=zeros(nx,ny);
amp520_d3=zeros(nx,ny);
amp660_d3=zeros(nx,ny);
d410_d3=zeros(nx,ny);
d520_d3=zeros(nx,ny);
d660_d3=zeros(nx,ny);
% 
for i=1:nx
    for j=1:ny
        d0=squeeze(d3d0_depth(:,i,j));
        keep1=find(depth0>=drange1(1) & depth0<=drange1(2));
        keep2=find(depth0>=drange2(1) & depth0<=drange2(2));
        keep3=find(depth0>=drange3(1) & depth0<=drange3(2));
        % find the maximum amplitude
        [amp410_d0(i,j),i1]=max(d0(keep1));
        [amp660_d0(i,j),i2]=max(d0(keep2));
        [amp520_d0(i,j),i3]=max(d0(keep3));
        ind410=keep1(i1);
        ind660=keep2(i2);
        ind520=keep3(i3);
        d410_d0(i,j)=depth0(ind410);
        d660_d0(i,j)=depth0(ind660);
        d520_d0(i,j)=depth0(ind520);

        d0=squeeze(d3d1_depth(:,i,j));
        % find the maximum amplitude
        [amp410_d1(i,j),i1]=max(d0(keep1));
        [amp660_d1(i,j),i2]=max(d0(keep2));
        [amp520_d1(i,j),i3]=max(d0(keep3));
        ind410=keep1(i1);
        ind660=keep2(i2);
        ind520=keep3(i3);
        d410_d1(i,j)=depth0(ind410);
        d660_d1(i,j)=depth0(ind660);
        d520_d1(i,j)=depth0(ind520);

        d0=squeeze(d3d2_depth(:,i,j));
        % find the maximum amplitude
        [amp410_d2(i,j),i1]=max(d0(keep1));
        [amp660_d2(i,j),i2]=max(d0(keep2));
        [amp520_d2(i,j),i3]=max(d0(keep3));
        ind410=keep1(i1);
        ind660=keep2(i2);
        ind520=keep3(i3);
        d410_d2(i,j)=depth0(ind410);
        d660_d2(i,j)=depth0(ind660);
        d520_d2(i,j)=depth0(ind520);

%         d0=squeeze(d3d3_depth(:,i,j));
%         % find the maximum amplitude
%         [amp410_d3(i,j),i1]=max(d0(keep1));
%         [amp660_d3(i,j),i2]=max(d0(keep2));
%         [amp520_d3(i,j),i3]=max(d0(keep3));
%         ind410=keep1(i1);
%         ind660=keep2(i2);
%         ind520=keep3(i3);
%         d410_d3(i,j)=depth0(ind410);
%         d660_d3(i,j)=depth0(ind660);
%         d520_d3(i,j)=depth0(ind520);
    end
end

d410_d0=d410_d0';
d520_d0=d520_d0';
d660_d0=d660_d0';

d410_d1=d410_d1';
d520_d1=d520_d1';
d660_d1=d660_d1';
% 
d410_d2=d410_d2';
d520_d2=d520_d2';
d660_d2=d660_d2';
% 
% d410_d3=d410_d3';
% d520_d3=d520_d3';
% d660_d3=d660_d3';
%% plot 410
dx=5; dy=5; dh=2;
xmin=110; ymin=20; hmin=100;
xmax=160; ymax=60; hmax=170;
% define grid center
x = xmin+dx/2:dx:xmax;
y = ymin+dy/2:dy:ymax;
h = hmin+dh/2:dh:hmax;

lonlim=[xmin xmax];
latlim=[ymin ymax];

[Xgrid,Ygrid]=meshgrid(x,y);

% save the model
% results=[Xgrid(:) Ygrid(:) d410_d0(:) d660_d0(:)];
% save Oboue22.dat results '-ascii'

Vgrid_d0=d410_d0;
F410_d0=scatteredInterpolant(Xgrid(:),Ygrid(:),Vgrid_d0(:),'natural','none');
Vgrid_d1=d410_d1;
F410_d1=scatteredInterpolant(Xgrid(:),Ygrid(:),Vgrid_d1(:),'natural','none');
Vgrid_d2=d410_d2;
F410_d2=scatteredInterpolant(Xgrid(:),Ygrid(:),Vgrid_d2(:),'natural','none');
% Vgrid_d3=d410_d3;
% F410_d3=scatteredInterpolant(Xgrid(:),Ygrid(:),Vgrid_d3(:),'natural','none');

xi=lonlim(1):0.2:lonlim(2);
yi=latlim(1):0.2:latlim(2);
[XI,YI]=meshgrid(xi,yi);
d410_d0_interp=F410_d0(XI,YI);
d410_d1_interp=F410_d1(XI,YI);
d410_d2_interp=F410_d2(XI,YI);
% d410_d3_interp=F410_d3(XI,YI);

% smooth the results
% ngrid=floor(5/0.2);
ngrid=1;
K = (1/ngrid^2)*ones(ngrid,ngrid);
d410_d0_smooth = conv2(d410_d0_interp,K,'same');
d410_d1_smooth = conv2(d410_d1_interp,K,'same');
d410_d2_smooth = conv2(d410_d2_interp,K,'same');
% d410_d3_smooth = conv2(d410_d3_interp,K,'same');

vmean=mean(d410_d1_smooth(~isnan(d410_d1_smooth)));

addpath 'm_map'

figure('units','normalized','Position',[0.0 0.0 1 1],'color','w'); 
set(gcf,'Position',[100 100 1000 1000],'color','w')
subplot(1,3,1)
m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
set(gcf,'color','w')
h=m_pcolor(XI,YI,d410_d0_smooth);
set(h,'edgecolor','none')
cm=colormap(jet);
colormap(flipud(cm));
caxis([385 435])
% 
m_gshhs('i','line','color','k','linewidth',1)
m_gshhs('lb2','line','color','k')
m_grid('linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',14);
title('MTZ 410 (True)','fontsize',18)
text(-0.16,0.98,'(a)','Units','normalized','FontSize',24)

subplot(1,3,2)
m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
set(gcf,'color','w')
h=m_pcolor(XI,YI,d410_d1_smooth);
set(h,'edgecolor','none')
cm=colormap(jet);
colormap(flipud(cm));
caxis([385 435])

m_gshhs('i','line','color','k','linewidth',1)
m_gshhs('lb2','line','color','k')
m_grid('linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',14);
title('MTZ 410 (Raw)','fontsize',18)
text(-0.1,0.98,'(b)','Units','normalized','FontSize',24)
% 
subplot(1,3,3)
m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
set(gcf,'color','w')
h=m_pcolor(XI,YI,d410_d2_smooth);
set(h,'edgecolor','none')
cm=colormap(jet);
colormap(flipud(cm));
caxis([385 435])

m_gshhs('i','line','color','k','linewidth',1)
m_gshhs('lb2','line','color','k')
m_grid('linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',14);
title('MTZ 410 (Reconstructed)','fontsize',18)
text(-0.16,0.98,'(c)','Units','normalized','FontSize',24)

% subplot(224)
% m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
% set(gcf,'color','w')
% h=m_pcolor(XI,YI,d410_d3_smooth);
% set(h,'edgecolor','none')
% cm=colormap(jet);
% colormap(flipud(cm));
% caxis([390 450])
% 
% m_gshhs('i','line','color','k','linewidth',1)
% m_gshhs('lb2','line','color','k')
% m_grid('linewidth',2,'tickdir','out',...
%     'xaxisloc','bottom','yaxisloc','left','fontsize',14);
% hh=colorbar('h');
% set(hh,'fontsize',14);
% set(hh,'Position',[0.3065 0.06 0.4220 0.0250])
% xlabel(hh,'Depth (km)');
% title('MTZ 410 depth (r=500,snr<=0.5)','fontsize',18)
% text(-0.16,0.98,'(d)','Units','normalized','FontSize',24)
% % export_fig 'mtz_410.png' '-r150'
% 
% % save to GMT plot
% outdir='/home/bm/Desktop/sspGMT';
% outdir='/Users/oboue/Desktop/ssp/synth_ssp_New/SS-recon_v0new/sspGMT';
% results=[Xgrid(:),Ygrid(:),Vgrid_d0(:)];
% save(fullfile(outdir,'d410true.txt'),'results','-ascii');
% 
% results=[Xgrid(:),Ygrid(:),Vgrid_d1(:)];
% save(fullfile(outdir,'d410raw.txt'),'results','-ascii');
% 
% results=[Xgrid(:),Ygrid(:),Vgrid_d2(:)];
% save(fullfile(outdir,'d410rec.txt'),'results','-ascii');
%% plot 520
% Vgrid_d0=d520_d0;
% F520_d0=scatteredInterpolant(Xgrid(:),Ygrid(:),Vgrid_d0(:),'natural','none');
% Vgrid_d1=d520_d1;
% F520_d1=scatteredInterpolant(Xgrid(:),Ygrid(:),Vgrid_d1(:),'natural','none');
% Vgrid_d2=d520_d2;
% F520_d2=scatteredInterpolant(Xgrid(:),Ygrid(:),Vgrid_d2(:),'natural','none');
% Vgrid_d3=d520_d3;
% F520_d3=scatteredInterpolant(Xgrid(:),Ygrid(:),Vgrid_d3(:),'natural','none');
% 
% xi=lonlim(1):0.2:lonlim(2);
% yi=latlim(1):0.2:latlim(2);
% [XI,YI]=meshgrid(xi,yi);
% d520_d0_interp=F520_d0(XI,YI);
% d520_d1_interp=F520_d1(XI,YI);
% d520_d2_interp=F520_d2(XI,YI);
% d520_d3_interp=F520_d3(XI,YI);
% 
% % ngrid=floor(2.5/0.2);
% ngrid=1;
% K = (1/ngrid^2)*ones(ngrid,ngrid);
% d520_d0_smooth = conv2(d520_d0_interp,K,'same');
% d520_d1_smooth = conv2(d520_d1_interp,K,'same');
% d520_d2_smooth = conv2(d520_d2_interp,K,'same');
% d520_d3_smooth = conv2(d520_d3_interp,K,'same');
% 
% vmean=mean(d520_d1_smooth(~isnan(d520_d1_smooth)));
% 
% figure;
% set(gcf,'Position',[100 100 1000 1000],'color','w')
% subplot(221)
% m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
% set(gcf,'color','w')
% h=m_pcolor(XI,YI,d520_d0_smooth);
% set(h,'edgecolor','none')
% cm=colormap(jet);
% colormap(flipud(cm));
% caxis([490 550])
% 
% m_gshhs('i','line','color','k','linewidth',1)
% m_gshhs('lb2','line','color','k')
% m_grid('linewidth',2,'tickdir','out',...
%     'xaxisloc','bottom','yaxisloc','left','fontsize',14);
% title('MTZ 520 depth (r=150,snr<=3)','fontsize',18)
% text(-0.16,0.98,'(a)','Units','normalized','FontSize',24)
% 
% subplot(222)
% m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
% set(gcf,'color','w')
% h=m_pcolor(XI,YI,d520_d1_smooth);
% set(h,'edgecolor','none')
% cm=colormap(jet);
% colormap(flipud(cm));
% caxis([490 550])
% 
% m_gshhs('i','line','color','k','linewidth',1)
% m_gshhs('lb2','line','color','k')
% m_grid('linewidth',2,'tickdir','out',...
%     'xaxisloc','bottom','yaxisloc','left','fontsize',14);
% title('MTZ 520 depth (r=150,snr<=0.5)','fontsize',18)
% text(-0.1,0.98,'(b)','Units','normalized','FontSize',24)
% 
% subplot(223)
% m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
% set(gcf,'color','w')
% h=m_pcolor(XI,YI,d520_d2_smooth);
% set(h,'edgecolor','none')
% cm=colormap(jet);
% colormap(flipud(cm));
% caxis([490 550])
% 
% m_gshhs('i','line','color','k','linewidth',1)
% m_gshhs('lb2','line','color','k')
% m_grid('linewidth',2,'tickdir','out',...
%     'xaxisloc','bottom','yaxisloc','left','fontsize',14);
% title('MTZ 520 depth (r=500,snr<=3)','fontsize',18)
% text(-0.16,0.98,'(c)','Units','normalized','FontSize',24)
% 
% subplot(224)
% m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
% set(gcf,'color','w')
% h=m_pcolor(XI,YI,d520_d3_smooth);
% set(h,'edgecolor','none')
% cm=colormap(jet);
% colormap(flipud(cm));
% caxis([490 550])
% 
% m_gshhs('i','line','color','k','linewidth',1)
% m_gshhs('lb2','line','color','k')
% m_grid('linewidth',2,'tickdir','out',...
%     'xaxisloc','bottom','yaxisloc','left','fontsize',14);
% hh=colorbar('h');
% set(hh,'fontsize',14);
% set(hh,'Position',[0.3065 0.06 0.4220 0.0250])
% xlabel(hh,'Depth (km)');
% title('MTZ 520 depth (r=500,snr<=0.5)','fontsize',18)
% text(-0.16,0.98,'(d)','Units','normalized','FontSize',24)
% export_fig 'mtz_520.png' '-r150'
% % save to GMT plot
% outdir='/Users/yunfeng/30_40/publications/oboue/ss/figures/mtz_map';
% results=[Xgrid(:),Ygrid(:),Vgrid_d1(:)];
% save(fullfile(outdir,'d520.txt'),'results','-ascii');
%% plot 660
Vgrid_d0=d660_d0;
F660_d0=scatteredInterpolant(Xgrid(:),Ygrid(:),Vgrid_d0(:),'natural','none');
Vgrid_d1=d660_d1;
F660_d1=scatteredInterpolant(Xgrid(:),Ygrid(:),Vgrid_d1(:),'natural','none');
Vgrid_d2=d660_d2;
F660_d2=scatteredInterpolant(Xgrid(:),Ygrid(:),Vgrid_d2(:),'natural','none');
% Vgrid_d3=d660_d3;
% F660_d3=scatteredInterpolant(Xgrid(:),Ygrid(:),Vgrid_d3(:),'natural','none');

xi=lonlim(1):0.2:lonlim(2);
yi=latlim(1):0.2:latlim(2);
[XI,YI]=meshgrid(xi,yi);
d660_d0_interp=F660_d0(XI,YI);
d660_d1_interp=F660_d1(XI,YI);
d660_d2_interp=F660_d2(XI,YI);
% d660_d3_interp=F660_d3(XI,YI);

K = (1/ngrid^2)*ones(ngrid,ngrid);
d660_d0_smooth = conv2(d660_d0_interp,K,'same');
d660_d1_smooth = conv2(d660_d1_interp,K,'same');
d660_d2_smooth = conv2(d660_d2_interp,K,'same');
% d660_d3_smooth = conv2(d660_d3_interp,K,'same');

vmean=mean(d660_d1_smooth(~isnan(d660_d1_smooth)));

figure('units','normalized','Position',[0.0 0.0 1 1],'color','w'); 
set(gcf,'Position',[100 100 1000 1000],'color','w')
subplot(1,3,1)
m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
set(gcf,'color','w')
h=m_pcolor(XI,YI,d660_d0_smooth);
set(h,'edgecolor','none')
cm=colormap(jet);
colormap(flipud(cm));
caxis([630 690])

m_gshhs('i','line','color','k','linewidth',1)
m_gshhs('lb2','line','color','k')
m_grid('linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',14);
title('MTZ 660 (True)','fontsize',18)
text(-0.16,0.98,'(a)','Units','normalized','FontSize',24)

subplot(1,3,2)
m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
set(gcf,'color','w')
h=m_pcolor(XI,YI,d660_d1_smooth);
set(h,'edgecolor','none')
cm=colormap(jet);
colormap(flipud(cm));
caxis([630 690])

m_gshhs('i','line','color','k','linewidth',1)
m_gshhs('lb2','line','color','k')
m_grid('linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',14);
title('MTZ 660 (Raw)','fontsize',18)
text(-0.1,0.98,'(b)','Units','normalized','FontSize',24)

subplot(1,3,3)
m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
set(gcf,'color','w')
h=m_pcolor(XI,YI,d660_d2_smooth);
set(h,'edgecolor','none')
cm=colormap(jet);
colormap(flipud(cm));
caxis([630 690])

m_gshhs('i','line','color','k','linewidth',1)
m_gshhs('lb2','line','color','k')
m_grid('linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',14);
 title('MTZ 660 (Reconstructed)','fontsize',18)
text(-0.16,0.98,'(c)','Units','normalized','FontSize',24)

% subplot(224)
% m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
% set(gcf,'color','w')
% h=m_pcolor(XI,YI,d660_d3_smooth);
% set(h,'edgecolor','none')
% cm=colormap(jet);
% cm=colormap(jet);
% colormap(flipud(cm));
% caxis([630 690])
% 
% m_gshhs('i','line','color','k','linewidth',1)
% m_gshhs('lb2','line','color','k')
% m_grid('linewidth',2,'tickdir','out',...
%     'xaxisloc','bottom','yaxisloc','left','fontsize',14);
% hh=colorbar('h');
% set(hh,'fontsize',14);
% set(hh,'Position',[0.3065 0.06 0.4220 0.0250])
% xlabel(hh,'Depth (km)');
% % title('MTZ 660 depth (r=500,snr<=0.5)','fontsize',18)
% text(-0.16,0.98,'(d)','Units','normalized','FontSize',24)
% export_fig 'mtz_660.png' '-r150'
% % save to GMT plot
% outdir='/Users/yunfeng/30_40/publications/oboue/ss/figures/mtz_map';
% results=[Xgrid(:),Ygrid(:),Vgrid_d1(:)];
% save(fullfile(outdir,'d660.txt'),'results','-ascii');


% % outdir='//home/bm/Desktop/sspGMT';
% outdir='/Users/oboue/Desktop/ssp/synth_ssp_New/SS-recon_v0new/sspGMT';
% results=[Xgrid(:),Ygrid(:),Vgrid_d0(:)];
% save(fullfile(outdir,'d660true.txt'),'results','-ascii');
% 
% results=[Xgrid(:),Ygrid(:),Vgrid_d1(:)];
% save(fullfile(outdir,'d660raw.txt'),'results','-ascii');
% 
% results=[Xgrid(:),Ygrid(:),Vgrid_d2(:)];
% save(fullfile(outdir,'d660rec.txt'),'results','-ascii');
%% plot thickness of MTZ (660-410)

thi_d0=d660_d0_smooth-d410_d0_smooth;
thi_d1=d660_d1_smooth-d410_d1_smooth;
thi_d2=d660_d2_smooth-d410_d2_smooth;
% thi_d3=d660_d3_smooth-d410_d3_smooth;
% 
ngrid=1;

K = (1/ngrid^2)*ones(ngrid,ngrid);
thi_d0_smooth = conv2(thi_d0,K,'same');
thi_d1_smooth = conv2(thi_d1,K,'same');
thi_d2_smooth = conv2(thi_d2,K,'same');
% thi_d3_smooth = conv2(thi_d3,K,'same');
vmean=mean(thi_d1_smooth(~isnan(thi_d1_smooth)));

figure('units','normalized','Position',[0.0 0.0 1 1],'color','w'); 
set(gcf,'Position',[100 100 1000 1000],'color','w')
subplot(1,3,1)
m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
set(gcf,'color','w')
h=m_pcolor(XI,YI,thi_d0_smooth);
set(h,'edgecolor','none')
cm=colormap(jet);
colormap(flipud(cm));
caxis([vmean-10 vmean+10])

m_gshhs('i','line','color','k','linewidth',1)
m_gshhs('lb2','line','color','k')
m_grid('linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',14);
title('Thickness (True)','fontsize',18)
text(-0.16,0.98,'(a)','Units','normalized','FontSize',24)

subplot(1,3,2)
m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
set(gcf,'color','w')
h=m_pcolor(XI,YI,thi_d1_smooth);
set(h,'edgecolor','none')
cm=colormap(jet);
colormap(flipud(cm));
caxis([vmean-10 vmean+10])

m_gshhs('i','line','color','k','linewidth',1)
m_gshhs('lb2','line','color','k')
m_grid('linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',14);
title('Thickness (Raw)','fontsize',18)
text(-0.1,0.98,'(b)','Units','normalized','FontSize',24)

subplot(1,3,3)
m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
set(gcf,'color','w')
h=m_pcolor(XI,YI,thi_d2_smooth);
set(h,'edgecolor','none')
cm=colormap(jet);
colormap(flipud(cm));
caxis([vmean-10 vmean+10])

m_gshhs('i','line','color','k','linewidth',1)
m_gshhs('lb2','line','color','k')
m_grid('linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',14);
title('Thickness (Reconstructed)','fontsize',18)
text(-0.16,0.98,'(c)','Units','normalized','FontSize',24)

% subplot(224)
% m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
% set(gcf,'color','w')
% h=m_pcolor(XI,YI,thi_d3_smooth);
% set(h,'edgecolor','none')
% cm=colormap(jet);
% colormap(flipud(cm));
% caxis([vmean-10 vmean+10])
% 
% m_gshhs('i','line','color','k','linewidth',1)
% m_gshhs('lb2','line','color','k')
% m_grid('linewidth',2,'tickdir','out',...
%     'xaxisloc','bottom','yaxisloc','left','fontsize',14);
% hh=colorbar('h');
% set(hh,'fontsize',14);
% set(hh,'Position',[0.3065 0.06 0.4220 0.0250])
% xlabel(hh,'Depth (km)');
% title('MTZ thickness (r=500,snr<=0.5)','fontsize',18)
% text(-0.16,0.98,'(d)','Units','normalized','FontSize',24)