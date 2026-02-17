function [aucc, pyrfracs, aucr, pyrfracsr, aucm]=P2L3sAUC(parms,fdv,rfdv)

%Unpack variables:
for ii=1:length(parms)
    x.(fdv.fitvars{ii})=parms(ii);
end
for ii=1:length(fdv.knowns)
    x.(fdv.knowns{ii})=fdv.knownvals(ii);
end

aucp=trapz(fdv.data(1,:).*fdv.TR);
aucl=trapz(fdv.data(2,:).*fdv.TR);
aucm=aucl/aucp;

vee=(1-x.vb)*x.vef;

% nbins=length(fdv.SPBins);
% pb=zeros(1,nbins);
% pe=zeros(1,nbins);
% pc=zeros(1,nbins);
% lc=zeros(1,nbins);

TR=fdv.TR(1); % assume constant TR
%Calculate flip angles for this part of slice profile:
FAp=fdv.FlipAngle(1,1)*fdv.SPBins;
FAl=fdv.FlipAngle(2,1)*fdv.SPBins;
%...and attenuation coefficients with averaged excitation losses:
apee=-((x.kve+x.kecp)/vee + 1/x.T1Pyr + (1-cosd(FAp))/TR);
apc=-(x.kpl + 1/x.T1Pyr + (1-cosd(FAp))/TR);
alc=-(1/x.T1Lac + (1-cosd(FAl))/TR + x.R1Lx);
%The AUC for each bin within the slice profile:
pb  =  sind(FAp)*exp(-fdv.TE/x.T2spv)*x.vb;
pe  = -sind(FAp)*exp(-fdv.TE/x.T2spe)*x.kve./apee;
pc  =  sind(FAp)*exp(-fdv.TE/x.T2spc)*x.kve*x.kecp/vee./apee./apc;
lc  = -sind(FAl)*exp(-fdv.TE/x.T2slc)*x.kve*x.kecp*x.kpl/vee./apee./apc./alc;

% for zz=1:length(fdv.SPBins)
%     %Calculate flip angles for this part of slice profile:
%     FAp=fdv.FlipAngle(1,1)*fdv.SPBins(zz);
%     FAl=fdv.FlipAngle(2,1)*fdv.SPBins(zz);
%     %...and attenuation coefficients with averaged excitation losses:
%     apee=-((x.kve+x.kecp)/vee + 1/x.T1Pyr + (1-cosd(FAp))/TR);
%     apc=-(x.kpl + 1/x.T1Pyr + (1-cosd(FAp))/TR);
%     alc=-(1/x.T1Lac + (1-cosd(FAl))/TR + x.R1Lx);
%     %The AUC for each bin within the slice profile:
%     pb(zz)  =  sind(FAp)*exp(-fdv.TE/x.T2spv)*x.vb;
%     pe(zz)  = -sind(FAp)*exp(-fdv.TE/x.T2spe)*x.kve/apee;
%     pc(zz)  =  sind(FAp)*exp(-fdv.TE/x.T2spc)*x.kve*x.kecp/vee/apee/apc;
%     lc(zz)  = -sind(FAl)*exp(-fdv.TE/x.T2slc)*x.kve*x.kecp*x.kpl/vee/apee/apc/alc;
% end
%figure(99);plot([pb;pe;pc;lc]');legend('pyr_v','pyr_e','pyr_c','lac_c');
aucc=sum(lc.*fdv.SPWeights)/sum((pb+pe+pc).*fdv.SPWeights);
pyrfracs = [sum(pb.*fdv.SPWeights) sum(pe.*fdv.SPWeights) sum(pc.*fdv.SPWeights)];
pyrfracs = pyrfracs/sum(pyrfracs);

if nargin>2
    % Recalculate with reference parameters:
    % nbins=length(rfdv.SPBins);
    % pb=zeros(1,nbins);
    % pe=zeros(1,nbins);
    % pc=zeros(1,nbins);
    % lc=zeros(1,nbins);
    
    TR=rfdv.TR(1);
    FAp=rfdv.FlipAngle(1,1)*rfdv.SPBins;
    FAl=rfdv.FlipAngle(2,1)*rfdv.SPBins;
    apee=-((x.kve+x.kecp)/vee + 1/x.T1Pyr + (1-cosd(FAp))/TR);
    apc=-(x.kpl + 1/x.T1Pyr + (1-cosd(FAp))/TR);
    alc=-(1/x.T1Lac + (1-cosd(FAl))/TR + x.R1Lx);
    pb  =  sind(FAp)*exp(-rfdv.TE/x.T2spv)*x.vb;
    pe  = -sind(FAp)*exp(-rfdv.TE/x.T2spe)*x.kve./apee;
    pc  =  sind(FAp)*exp(-rfdv.TE/x.T2spc)*x.kve*x.kecp/vee./apee./apc;
    lc  = -sind(FAl)*exp(-rfdv.TE/x.T2slc)*x.kve*x.kecp*x.kpl/vee./apee./apc./alc;
    % for zz=1:length(rfdv.SPBins)
    %     FAp=rfdv.FlipAngle(1,1)*rfdv.SPBins(zz);
    %     FAl=rfdv.FlipAngle(2,1)*rfdv.SPBins(zz);
    %     apee=-((x.kve+x.kecp)/vee + 1/x.T1Pyr + (1-cosd(FAp))/TR);
    %     apc=-(x.kpl + 1/x.T1Pyr + (1-cosd(FAp))/TR);
    %     alc=-(1/x.T1Lac + (1-cosd(FAl))/TR + x.R1Lx);
    %     pb(zz)  =  sind(FAp)*exp(-rfdv.TE/x.T2spv)*x.vb;
    %     pe(zz)  = -sind(FAp)*exp(-rfdv.TE/x.T2spe)*x.kve/apee;
    %     pc(zz)  =  sind(FAp)*exp(-rfdv.TE/x.T2spc)*x.kve*x.kecp/vee/apee/apc;
    %     lc(zz)  = -sind(FAl)*exp(-rfdv.TE/x.T2slc)*x.kve*x.kecp*x.kpl/vee/apee/apc/alc;
    % end
    aucr=sum(lc.*rfdv.SPWeights)/sum((pb+pe+pc).*rfdv.SPWeights);
    %pyrfracs=[pfbc pfec pfcc]/(pfbc+pfec+pfcc);
    pyrfracsr = [sum(pb.*rfdv.SPWeights) sum(pe.*rfdv.SPWeights) sum(pc.*rfdv.SPWeights)];
    pyrfracsr = pyrfracsr/sum(pyrfracsr);
else
    aucr=[];
    pyrfracsr=[];
end