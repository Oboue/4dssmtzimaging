function [d,is_reversal] = ss_check_polarity(d,t,tss)
% check the ampltidue of the SS phase and make sure that the absolute value
% the maximum possitive amplitude is always larger than that of the negative
% amplitude, i.e., abs(max(sig)) > abs(min(sig)).
% This follows the criteria proposed by Peter Shearer?
% Input:
% d: data
% t: time
% tss: arrival time of SS phase
% swin: signal window length
% Output:
% d: data
% is_flip: flag that is true for polarity reversal and false otherwise
% default value
if nargin<4
    swin=100;
end
is_reversal=false;
% index of SS main phase
dt=t(2)-t(1);
ii=floor((tss-t(1)+dt/2)/dt);
is=floor(ii-(swin/dt)/2:ii+(swin/dt)/2);
sig=d(is);
if abs(max(sig)) < abs(min(sig))
    d = -d;
    is_reversal=true;
end