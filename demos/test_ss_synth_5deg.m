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
addpath('../etopo1_bed_c_f4') 
javaaddpath('./FMI/lib/FMI.jar');
addpath('~/MATLAB/m_map');
addpath('../Meier_2009') 
addpath('../plot_tectonic_fold_maps') 
addpath('../slab')
addpath('../MatSAC');
addpath('../utils');
addpath('../subroutines');
addpath('./FMI/matTaup');
addpath('../TX2019slab')
addpath('../irisFetch-matlab-2.0.12')
addpath export_fig
% addpath /Users/oboue/Desktop/desktop03-02-2023/ssp/ssreconsynth
%% load the data
% datadir = './DATA/s20_topo/';
% sacfiles = dir(fullfile(datadir,'*sac'));
% ss=[];
% for i = 1:length(sacfiles) 
%     if mod(i,100) == 0
%         disp(['Reading ',num2str(i),'/',num2str(length(sacfiles)),' traces']);
%     end
%     [t,d,SAChdr] = fget_sac(fullfile(sacfiles(i).folder,sacfiles(i).name));
% %     tmp.d = d;
% %     tmp.t = t;
%     tmp.stla=SAChdr.station.stla;
%     tmp.stlo=SAChdr.station.stlo;
%     tmp.stel=SAChdr.station.stel;
%     tmp.sta =SAChdr.station.kstnm;
%     tmp.evla=SAChdr.event.evla;
%     tmp.evlo=SAChdr.event.evlo;
%     tmp.evdp=SAChdr.event.evdp/1000.; % meter to kilometer
%     tmp.mag=SAChdr.event.mag;
%     tmp.dist=SAChdr.evsta.dist;
%     tmp.az=SAChdr.evsta.az;
%     tmp.baz=SAChdr.evsta.baz;
%     tmp.gcarc=SAChdr.evsta.gcarc;
%     % calculate the bounce point (midpoint)
%     [tmp.bplat,tmp.bplon]=gc_midpoint(tmp.evla, tmp.evlo, tmp.stla, tmp.stlo);
%     % calculate SS arrival time
%     times=taupTime('prem',tmp.evdp,'SS,S^410S,S^660S','sta',[tmp.stla tmp.stlo],...
%         'evt',[tmp.evla,tmp.evlo]);
%     tmp.t660=times(1).time;
%     tmp.t410=times(2).time;
%     tmp.tss=times(3).time;
%     tmp.gcarc=times(1).distance;
%     t1=tmp.tss-900;
%     t2=tmp.tss+300;
%     [d_cut,t_cut]=chopSeis(d,t,t1,t2);
%     [d_resampled,dt,t_resampled]=resampleSeis(d_cut,t_cut,1);
% %     plot(t,d); hold on; plot(t_resampled,d_resampled)
%     tmp.d=d_resampled;
%     tmp.t=t_resampled;
%     ss=[ss tmp];
% end
% t=t_resampled;
% % save as mat
% save 'ss.mat' 'ss';

load ss_synth.mat;
%% calculate SNR
for k=1:length(ss)
    d=ss(k).d;
    t=ss(k).t;
    tss=ss(k).tss;
    ss(k).snr = ss_snr(d,t,tss);
end
%% check polarity reversal
for k=1:length(ss)
    d=ss(k).d;
    t=ss(k).t;
    tss=ss(k).tss;
    [d,is_reversal] = ss_check_polarity(d,t,tss);
    ss(k).is_reversal = is_reversal;
    if is_reversal
        ss(k).d = d;
    end
end
remove = [ss.snr]<=0;
ss(remove) = [];
%% apply cross-correlation to all traces
disp('Align traces')
nt=length(t);
for k=1:length(ss)
    ss(k).d=ss(k).d(1:nt);
end
din = [ss.d]; % flatten the tensor
N=5; % number of iteration for cross-correlation measurments
t=0:nt-1;
times = repmat(t(:),1,size(din,2));
t0=ones(1,size(din,2))*900; % SS arrival
xwin=[-100 100]; % cross-correlation window
maxlag=50; % maximum time lag
is_plot=1; % flag controls plotting
dout = ss_align_v2(din,times,N,t0,xwin,maxlag,is_plot);
for k=1:length(ss)
    ss(k).d=dout(:,k);
end
%% plot a few traces
% figure;
% for n=1:length(ss)
%     plot(ss(n).t,ss(n).d/max(ss(n).d)); hold on;
%     plot([ss(n).t410 ss(n).t410],[-1 1],'r--')
%     plot([ss(n).t660 ss(n).t660],[-1 1],'r--')
%     pause;
%     clf
% end
%% binning
disp('Binning')
% dx=2.5; dy=2.5; dh=2;
dx=5; dy=5; dh=2;

xmin=110; ymin=20; hmin=100;
xmax=160; ymax=60; hmax=170;
% define grid center
x = xmin+dx/2:dx:xmax;
y = ymin+dy/2:dy:ymax;
h = hmin+dh/2:dh:hmax;
t = ss(1).t;
lonlim=[xmin xmax];
latlim=[ymin ymax];
nx=length(x); ny=length(y); nh=length(h); nt = length(t);
% 5D case
% d1=zeros(nt, nx, ny, nh, nphi);
% 4D case
d1 = zeros(nt, nx, ny, nh);
fold_map=zeros(nx,ny,nh);
flow=1/75.;
fhigh=1/15.;
for n=1:length(ss)
    j=floor((ss(n).bplat-ymin)/dy)+1;
    i=floor((ss(n).bplon-xmin)/dx)+1;
    k=floor((ss(n).gcarc-hmin)/dh)+1;
    if i<=0 || i>nx || j<=0 || j>ny || k<=0 || k>nh
        continue;
    end
%     l=floor(ss(n).phi/dphi)+1;
    fold_map(i,j,k)=fold_map(i,j,k)+1;
    d=ss(n).d;
    % bandpass filter
    d_filt=bandpassSeis(d,1,flow,fhigh,3);
    d_filt=d_filt/max(d_filt);
    d1(:,i,j,k)=d1(:,i,j,k)+d_filt(1:nt);
end
% nomalization
for i=1:nx
    for j=1:ny
        for k=1:nh
            if fold_map(i,j,k)>0
               d1(:,i,j,k)=d1(:,i,j,k)/fold_map(i,j,k); 
            end
        end
    end
end

% compute the theoretical precursor arrival times
for n=1:nh
    times=taupTime('ak135',10,'SS,S^410S,S^660S','deg',h(n));
    indices = find(strcmp({times.phaseName},'S^660S'));
    t660(n)=times(indices(1)).time;
    indices = find(strcmp({times.phaseName},'S^410S'));
    t410(n)=times(indices(1)).time;
    indices = find(strcmp({times.phaseName},'SS'));
    tss(n)=times(indices(1)).time;
end

% stack to form a 2D section
d2d = squeeze(mean(mean(d1,3),2)); % simple averaging
W = any(d1);              % obtain the non-zero trace
w = squeeze(sum(sum(W,3),2));  % calcualte the weight
d2d_w = squeeze(sum(sum(d1,3),2))*diag(1./w); % weighted averaging

% find the time of SS phase and set it to 0 time
[~,index] = max(sum(d2d,2));
tshift = t(index);
t=t-tshift;
ntraces = squeeze(sum(sum(fold_map,2),1));
% plot 2D data
% figure;
% subplot(511)
% bar(h,ntraces)
% subplot(5,1,2:5)
% set(gcf,'Position',[0 0 1000 1000],'Color','w')
% wigb(d2d,10,h,t); hold on;
% plot(h,t660-tss,'--r')
% plot(h,t410-tss,'--r')
% axis xy
% ylim([-500 100])
% ylabel('Time (s)')
% xlabel('Distance (deg)')
% set(gca,'fontsize',14)
%% apply NMO correction to all traces
disp('Normal move-out correction')
d1_nmo=zeros(size(d1));
h0=135;
is_plot=false;
% parfor i=1:nx
for i=1:nx
    for j=1:ny
        for k=1:nh
            din = d1(:,i,j,k);
            if any(d)
                [dout,t410_ref,t660_ref] = ss_nmo(din,t,h(k),h0,is_plot);
                d1_nmo(:,i,j,k)=dout;
            end
        end
    end
end

% perform stacking for each CMP gather
d3d_nmo = zeros(nt,nx,ny);
d3d = zeros(nt,nx,ny);
for i=1:nx
    for j=1:ny
        % move-out corrected cmp
        d_cmp = squeeze(d1_nmo(:,i,j,:));
        nstack = sum(any(d_cmp));
        if nstack>0
            d_stack = sum(d_cmp,2)/nstack;
            d3d_nmo(:,i,j)=d_stack;
        end
        % non-move-out corrected cmp
        d_cmp = squeeze(d1(:,i,j,:));
        nstack = sum(any(d_cmp));
        if nstack>0
            d_stack = sum(d_cmp,2)/nstack;
            d3d(:,i,j)=d_stack;
        end
    end
end
% compare move-out corrected and non-move-out corrected CMP gather
% figure
% subplot(211)
% imagesc(1:nx*ny,t,reshape(d3d,nt,nx*ny))
% colormap(seismic(3))
% ylim([-500 100])
% axis xy;
% caxis([-0.05,0.05])
% subplot(212)
% imagesc(1:nx*ny,t,reshape(d3d_nmo,nt,nx*ny))
% colormap(seismic(3))
% ylim([-500 100])
% axis xy;
% caxis([-0.05,0.05])
%% 3D post-stack reconstruction using RDRR algorithm
disp('3D RDRR')
t0=-300;
t1=-100;
%t0=-500;
%t1=100;
keep = t>=t0 & t<=t1;
d3d=d3d_nmo(keep,:,:);
t=t(keep);
nt=length(t);

mask = repmat(any(d3d),size(d3d,1),1);

d3d=d3d.*mask;
dt=t(2)-t(1);
flow=1/75;
fhigh=1/15.;
Niter=10;
mode=1;
verb=1;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing

N=16;      % Rank 
K=3;       % Damping factor
u=0.00001; % Cooling factor 
e=0.9;     % Rational transfer function coefficient 
ws=1;      % Window size

% d3d_recon=fxyrdrr3d_denoising_recon(d3d,mask,flow,fhigh,dt,N,K,Niter,eps,verb,mode,a,u,e,ws);
% save output/d3reconfield.mat d3d d3d_recon x y t

% %% plot 3D reconstruction
% % figure;
% % set(gcf,'Position',[100 100 1600 500],'color','w')
% % imagesc(1:nx*ny,t,reshape(d3d,nt,nx*ny)); hold on;
% % plot(1:nx*ny,t410_ref*ones(1,nx*ny),'--r','linewidth',4)
% % plot(1:nx*ny,t660_ref*ones(1,nx*ny),'--r','linewidth',4)
% % title('Post-stack (raw)')
% % colormap(seis);caxis([-0.08 0.08]); colorbar;
% % axis xy
% % ylim([-300 -100])
% % ylabel('Time to SS (sec)')
% % xlabel('Trace');
% % set(gca,'fontsize',30)
% % 
% % figure;
% % set(gcf,'Position',[100 100 1600 500],'color','w')
% % imagesc(1:nx*ny,t,reshape(d3d_recon,nt,nx*ny)); hold on;
% % plot(1:nx*ny,t410_ref*ones(1,nx*ny),'--r','linewidth',4)
% % plot(1:nx*ny,t660_ref*ones(1,nx*ny),'--r','linewidth',4)
% % title('Post-stack (Reconstructed)')
% % colormap(seis);caxis([-0.08 0.08]); colorbar;
% % axis xy
% % ylim([-300 -100])
% % ylabel('Time to SS (sec)')
% % xlabel('Trace');
% % set(gca,'fontsize',30)
%% 4D reconstruction using RDRR algorithm
disp('4D RDRR')
d4d = d1_nmo(keep,:,:,:);
[nt,nx,ny,nh]=size(d4d);

mask = repmat(any(d4d),size(d4d,1),1);%
d4d=d4d.*mask;

% N=250; % Rank 
% K=1.5;  % Damping factor default 2 larger value preserve more signal 
% u=0.001; % Cooling factor level of sparsity can be fixed
% e=0.80;  % Rational transfer function coefficient can be fixed
% ws=1;    % Window size 1 by default
% iflb=0;
% d4d_recon=rdrr5d_lb_recon(d4d,mask,flow,fhigh,dt,N,K,Niter,eps,verb,mode,iflb,a,u,e,ws);
% % save d4d_recon.mat d4d_recon

% N=30; % Rank 
% K=2;  % Damping factor default 2 larger value preserve more signal 
% u=0.01; % Cooling factor level of sparsity can be fixed
% e=1;  % Rational transfer function coefficient can be fixed
% ws=1;    % Window size 1 by default
% iflb=0;
% d=rdrr5d_lb_recon(d4d_recon,mask,flow,fhigh,dt,N,K,Niter,eps,verb,mode,iflb,a,u,e,ws);
% % save synthssp4d.mat d
%% 4D synthetic ssp reconstruction
load synthssp4d.mat

% figure;
% subplot(1,1,1);
% set(gcf,'Position',[100 100 1600 500],'color','w')
% imagesc(1:nx*ny*nh,t,reshape(d,nt,nx*ny*nh)); hold on;
% plot(1:nx*ny*nh,t410_ref*ones(1,nx*ny*nh),'--r','linewidth',4)
% plot(1:nx*ny*nh,t660_ref*ones(1,nx*ny*nh),'--r','linewidth',4)
% title('True')
% colormap(seis);caxis([-0.12 0.12]); colorbar;
% axis xy
% ylim([-300 -100])
% ylabel('Time to SS (sec)')
% xlabel('Trace');
% set(gca,'fontsize',30)

% add random noise
randn('state',201314);

var=0.1;
synthssp4dn=d+var*randn(size(d));

ratio=0.3; % 70 percent missing traces 
mask=genmask(reshape(synthssp4dn,nt,nx*ny*nh),ratio,'c',201415);
mask=reshape(mask,nt,nx,ny,nh);
% size(ss4din)
synthssp4d0=synthssp4dn.*mask;
% save synthssp4d0.mat synthssp4d0 t410_ref t660_ref nt nx ny nh t mask
% o=yc_snr(d,synthssp4dn,2)
% load synthssp4d0
% o=yc_snr(d,synthssp4d0,2)

load synthssp4d.mat

load synthssp4d0.mat 

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