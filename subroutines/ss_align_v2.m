function [dout,tlag] = ss_align_v2(din,times,N,t0,xwin,maxlag,is_plot)
% cross-correlate traces to obtain the optimal time lag between the
% reference trace and each individual trace
% Input:
% din: nt*nx matrix that contains all traces
% times: nt*nx matrix that contains time axis
% N: number of iteration for cross-correlation measurements
% t0: nx*1 vector estimated arrival time of targeting phases
% xwin: 2*1 vector that defines the time window for xcorr
% maxlag: maximum time lag for xcorr
% is_plot: flag that controls plotting
% Output:
% dout: nt*nx matrix that contains timeshifted traces
% tlag: time difference between each trance and the reference trace
% change log: Jan. 1, 2022, output the time lag
if nargin < 7
    is_plot = 0;
end
trace_ref=[];
tlag=[];
% conduct N cross-correlation
for ncorr = 1:N
    keep = any(din);
    traces = din(:,keep)*diag(1./max(din(:,keep)));
    % cut the waveform arround P
    ntr = size(traces,2);
    tp=zeros(ntr,ntr);
    traces_cut=[];
    times_cut=[];
    for n = 1:ntr
        start = t0(n)+xwin(1);
        finish = t0(n)+xwin(2);
        [traces_cut(:,n), times_cut(:,n)] = chopSeis( traces(:,n), times(:,n), start, finish);
        traces_cut(:,n) = traces_cut(:,n)/max(traces_cut(:,n));
        % time difference in t0
        for m = 1:ntr
            tp(n,m)=t0(m)-t0(n);
        end
    end
    % reference trace
    trace_ref(:,ncorr) = mean(traces_cut,2);
    dt = times(2,1)-times(1,1);
    nlag = round(maxlag/dt);
    for n=1:ntr
        [C,lags] = xcorr(trace_ref(:,ncorr),traces_cut(:,n),nlag);
        [~,index] = max(C);
        % find the optimal delay time
        tdiff(n) = lags(index)*dt;
    end
    tlag=[tlag; tdiff];
    % apply timeshift to the traces
    traces=din(:,keep);
    traces_shift=[];
    for n = 1:length(tdiff)
        data = traces(:,n);
        time = times(:,n);
        datashift = fftShift(data,time,tdiff(n));
        traces_shift(:,n) = datashift;
    end
    if is_plot
        figure;
        imagesc(1:size(traces,2),time,[traces traces_shift])
        colormap(seismic(3))
        caxis([-0.0001,0.0001])
%         caxis([-0.05,0.05])
    end
    % update din
    din(:,keep)=traces_shift;
end
dout = din;
tlag=sum(tlag);