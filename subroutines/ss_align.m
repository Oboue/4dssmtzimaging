function dout = ss_align(din,times,N,t0,xwin,maxlag,is_plot)
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
if nargin < 7
    is_plot = 0;
end
trace_ref=[];
traces=[];
% conduct N cross-correlation
for ncorr = 1:N
    keep = any(din);
    din(:,keep) = din(:,keep)*diag(1./max(din(:,keep)));
    traces=din(:,keep);
    % reference trace
    trace_ref(:,ncorr) = mean(traces,2);
    % insert the reference trace
    traces = [trace_ref(:,ncorr) traces];
    
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
    
    dt = times(2,1)-times(1,1);
    maxlag = round(maxlag/dt);
    [C,lags] = xcorr(traces_cut,maxlag);
    [~,index] = max(C);
    I = reshape(index,ntr,ntr)';
    % find the optimal delay time
    tdelay = (lags(I(:)))*dt;
    tdelay = triu(reshape(tdelay,ntr,ntr)-tp);
    % take the first row, which is the time delay relative to the reference
    % trace
    tdiff = tdelay(1,:);
    
    % apply timeshift to the traces
    traces_shift=[];
    for n = 1:length(tdiff)
        data = traces(:,n);
        time = times(:,n);
        datashift = fftShift(data,time,tdiff(n));
        traces_shift(:,n) = datashift;
    end
    % remove the refernce trace
    traces(:,1)=[];
    traces_shift(:,1)=[];
    if is_plot
        figure;
        imagesc(1:size(traces,2),time,[traces traces_shift])
        colormap(seismic(3))
        caxis([-0.05,0.05])
    end
    % update din
    din(:,keep)=traces_shift;
end
dout = din;