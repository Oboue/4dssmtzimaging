function [bootdata,bootmean,bootsigma]=ss_bootstrap(d,nboot,ntrace)
% apply bootstrap to each cmp gather
% Input:
% d: data
% nboot: number of bootstrap iterations
% ntrace: number of traces that drawn from the data
% Output:
% bootdata: bootstrap data of size nt*nboot
% bootmean: mean value of bootstrap data
% bootsigma: standard deviation of bootstrap data
[nt,nx]=size(d);
bootdata=zeros(nt,nboot);
if nx == 1
    bootmean = d(:);
    bootsigma = zeros(nt,1);
    for n = 1:nboot
        bootdata(:,n) = bootmean;
    end
else
    % conduct resampling for nboot times
    for k = 1:nboot
        % everytime draw n traces randomally from the data (repeated samples are allowed)
        rng('shuffle','twister')
        indices = ceil(nx.*rand(ntrace,1));
        % calculate the mean (stack) of the picked traces and store the stack
        bootdata(:,k) = ignoreNaN(d(:,indices),@mean,2);
    end
    % calcualte the mean of stacks from all iterations (a stack of stacks)
    % and the corresponding standard deviation
    bootmean = ignoreNaN(bootdata',@mean);
    bootsigma = ignoreNaN(bootdata',@std);
end