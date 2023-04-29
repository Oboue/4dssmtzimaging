function snr = ss_snr(d,t,tss,swin,nwin)
% calculate the snr of SS phase, which is defined by the ratio of the maximum
% amplitude of the signal window to that of the noise window
% Input:
% d: data
% t: time
% tss: arrival time of SS phase
% swin: signal window length
% nwin: noise window length
% Output:
% snr: signal to noise ratio
% default value
if nargin<4
    swin=60;
    nwin=400;
end
dt=t(2)-t(1);
% index of SS main phase
ii=floor((tss-t(1)+dt/2)/dt);
is=floor(ii-(swin/dt)/2:ii+(swin/dt)/2);
in=floor(is(1)-1-nwin/dt:is(1)-1);
sig=d(is);
noi=d(in);
snr=max(abs(sig))/max(abs(noi));
% figure;
% plot(t,d); hold on; plot(t(is),d(is),'r'); plot(t(in),d(in),'b');