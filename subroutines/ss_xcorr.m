function dout = ss_xcorr(d,t,tss,swin)
% cross correlation to remove the source effect
% Input:
% d: data
% t: time
% tss: arrival time of SS phase
% swin: signal window length
% Output:
% snr: signal to noise ratio
% default value
if nargin<4
    swin=100;
end
dt=t(2)-t(1);
nt=length(t);
% index of SS main phase
ii=floor((tss-t(1)+dt/2)/dt);
is=floor(ii-(swin/dt)/2:ii+(swin/dt)/2);

sig=d(is);
% figure;
% plot(t,d); hold on; plot(t(is),d(is),'r'); plot(t(in),d(in),'b');
sig=taper(sig,20,20);
[xc,lag] = xcorr(d,sig);
imid=find(lag == 0);
% [~,imax]=max(xc);
i1=floor((tss-t(1))/dt);
i2=nt-i1-1;
imax=imid+i1-floor((swin/2/dt));
% tout=lag(imax-900:imax+400);
dout=xc(imax-i1:imax+i2);
dout=dout/max(dout);
[xc,lag]=xcorr(d,dout);
[~,index] = max(xc);
tdiff = lag(index)*dt;
dout = fftShift(dout,t,tdiff);

% figure;
% plot(d); hold on; 
% plot(dout);

