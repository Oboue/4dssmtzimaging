% Jun. 11, 2022, Yunfeng Chen, Global Seismology Group, Zhejiang University
% heterogeneity correction of ss precursor
function [z,tc,tt1d,tt3d] = ss_mancor(evla, evlo, stla, stlo, mod3d)
if nargin<5
    mod3d='S20RTS';
end
dcdps=[0,100:25:1000];
tSdS1d=zeros(size(dcdps));
tSdS3d=zeros(size(dcdps));
tcorr=zeros(size(dcdps));
for n=1:length(dcdps)
    dcdp=dcdps(n);
    disp(['Calculating depth ',num2str(dcdp)])
    dc=prem('depths',dcdp);
    dcvp=dc.vp;
    dcvs=dc.vs;
    dcrho=dc.rho;
    dcqk=dc.qk;
    dcqu=dc.qu;

    model=prem();
    depth=model.depth;
    % check if discontinuity exists in the model
    if any(dcdp==depth)
        dcdp=dcdp+0.1;
    end
    phstr=['S^',num2str(dcdp),'S'];
    taupobjfile=[pwd,'/taup_models/prem_dc',num2str(dcdp),'.taup'];
    if exist(taupobjfile)
        paths=tauppath('mod',taupobjfile,'dep',0,'ph',['SS,',phstr], ...
            'sta',[stla stlo], 'evt',[evla evlo]);
    else
        % insert the discontinuity into the model
        indx=find(dcdp>depth,1,'last');
        vp=model.vp;
        vs=model.vs;
        rho=model.rho;
        qk=model.qk;
        qu=model.qu;

        depth1=[depth(1:indx,1);dcdp;dcdp;depth(indx+1:end,1)];
        vp1=[vp(1:indx,1);dcvp;dcvp*(1.00001);vp(indx+1:end,1)];
        vs1=[vs(1:indx,1);dcvs;dcvs*(1.00001);vs(indx+1:end,1)];
        rho1=[rho(1:indx,1);dcrho;dcrho*(1.00001);rho(indx+1:end,1)];
        qk1=[qk(1:indx,1);dcqk;dcqk;qk(indx+1:end,1)];
        qu1=[qu(1:indx,1);dcqu;dcqu;qu(indx+1:end,1)];
        model.depth=depth1;
        model.vp=vp1;
        model.vs=vs1;
        model.rho=rho1;
        model.qk=qk1;
        model.qu=qu1;
        paths=tauppath('mod',model,'dep',0,'ph',['SS,',phstr], ...
            'sta',[stla stlo], 'evt',[evla evlo]);
        % write_1dmodel_nd('prem_dc.nd',model);
    end
    % do not consider the crust for mantle correction
    remove=[paths.distance]>180;
    paths(remove)=[];
    paths_new=crustless_raypaths(paths,'crust1.0');
    [corr,tt1d,tt3d]=mancor(paths_new,mod3d);
    % T_(SS-SdS)=(T^PREM_SS-T^TOMO_SS)-(T^PREM_SdS-T^TOMO_SdS)
    tSdS1d(n)=tt1d(2)-tt1d(1);
    tSdS3d(n)=tt3d(2)-tt3d(1);
    tcorr(n)=(tt1d(2)-tt3d(2))-(tt1d(1)-tt3d(1));
end
% interpolate to get all depths
z=0:1:1000;
tc=interp1(dcdps,tcorr,z);
tt1d=interp1(dcdps,tSdS1d,z);
tt3d=interp1(dcdps,tSdS3d,z);
% plot(tcorr,dcdps,'.'); hold on;
% plot(tc,z)