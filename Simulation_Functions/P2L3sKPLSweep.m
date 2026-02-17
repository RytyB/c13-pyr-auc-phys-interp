function [kplo, parmso, resido, residvec, kploi, nparmtrends]=P2L3sKPLSweep(fdv3s,kplvec,nguessperkpl,LB,UB,LGB,UGB)

verbose=1;

nkpls=length(kplvec);
nparms=length(fdv3s.fitvars);
kpli=find(strcmp('kpl',fdv3s.knowns));
kvei=find(strcmp('kve',fdv3s.fitvars));
vbi=find(strcmp('vb',fdv3s.fitvars));
vefi=find(strcmp('vef',fdv3s.fitvars));
kecpi=find(strcmp('kecp',fdv3s.fitvars));
R1Lxi=find(strcmp('R1Lx',fdv3s.fitvars));
VIFScalei=find(strcmp('VIFScale',fdv3s.fitvars));
vbi=find(strcmp('vb',fdv3s.fitvars));
linscale=1.0-strncmp(fdv3s.fitvars,'k',1)*1.0;

residvec=zeros(1,nkpls);
kplcsfits=zeros(nkpls,nparms);
kves=zeros(1,nkpls);
vbs=zeros(1,nkpls);
vefs=zeros(1,nkpls);
kecps=zeros(1,nkpls);
R1Lxs=zeros(1,nkpls);

for kplvi=1:nkpls 
    fprintf('%s - P2L3s - %s - kpl %d/%d\n',datetime('now'),fdv3s.Name,kplvi,nkpls);
    fdv3s.knownvals(kpli)=kplvec(kplvi);
    bestresids=ones(1,nguessperkpl)*Inf;
    bestfits=zeros(nguessperkpl,nparms);
    parfor fiti=1:nguessperkpl 
        linguess=LGB+rand(1,nparms).*(UGB-LGB);
        expguess=10.^(log10(LGB)+rand(1,nparms).*(log10(UGB)-log10(LGB)));
        expguess(isnan(expguess))=0; % in case the guess is zero and log(0)=NaN
        Guess=linguess.*linscale+expguess.*(1-linscale);
        %Guess(isnan(Guess))=0;
        [pm,pmi]=max(fdv3s.data(1,:));
        %[pm Guess(vbi) max(fdv3s.VIFP) fdv3s.FlipAngle(1,pmi)]
        Guess(VIFScalei)=(1+rand(1))*pm/Guess(vbi)/max(fdv3s.VIFP)/sind(fdv3s.FlipAngle(1,pmi));
        try
            [fits,resid] = lsqnonlin(@(x) P2L3sErr(x,fdv3s),Guess,LB,UB,optimset('display','off'));
        catch materr
            fprintf('Warning: %s\n',materr.identifier)
            fits=Guess;
            resid=sum(P2L3sErr(Guess,fdv3s).^2);
            if isnan(resid)
                resid=Inf;
            end
        end
        if resid<bestresids(fiti)
            bestfits(fiti,:)=fits;
            bestresids(fiti)=resid;
        end
    end % nguessperkpl
    [bestresid, bestresidi]=min(bestresids);
    residvec(kplvi)=bestresid;
    kplcsfits(kplvi,:)=bestfits(bestresidi,:);                
    vbs(kplvi)=bestfits(bestresidi,vbi);                   
    vefs(kplvi)=bestfits(bestresidi,vefi);         
    kves(kplvi)=bestfits(bestresidi,kvei);                           
    kecps(kplvi)=bestfits(bestresidi,kecpi);
    % vefs(kplvi)=bestfits(bestresidi,vefi); 
    R1Lxs(kplvi)=bestfits(bestresidi,R1Lxi);
    fprintf('P2L3s kpl = %8.3e --- bestresid = %7.5f\n\n',kplvec(kplvi),bestresid)
end % kplvec
%first find kpl value that is global min or 1% above:   
[kpl, parms, resid]=estparms(residvec,kplvec,kplcsfits,fdv3s,0);
[minresid,minresidi]=min(residvec);

nparmtrends=[vbs;kves;vefs;kecps;R1Lxs];

%don't search below kpl-1%, where we clearly dont have a solution:
kplli=find(kplvec>kpl,1); %kpl lower limit index

sosnp=sqrt(kves.*kves+kecps.*kecps+R1Lxs.*R1Lxs+(vefs-0.5).*(vefs-0.5));
[minsos minsosi]=min(sosnp(kplli:end));
minsosi=minsosi+kplli-1;

nkves=kves/UB(kvei);
nkecps=kecps/UB(kecpi);
nR1Lxs=R1Lxs/UB(R1Lxi);
nsosnp=sqrt(nkves.*nkves+nkecps.*nkecps+nR1Lxs.*nR1Lxs+(vefs-0.5).*(vefs-0.5));
[minnsos minnsosi]=min(nsosnp(kplli:end));
minnsosi=minnsosi+kplli-1;

%Note kplvec here is not normalized... assume values over ~1 are improbable
nsosnp2=sqrt((kplvec.^2)+(nkves.^2)+(nkecps.^2)+(nR1Lxs.^2)+((vefs-0.5).^2));
[minnsos2 minnsosi2]=min(nsosnp2(kplli:end));
minnsosi2=minnsosi2+kplli-1;

kplo=[kpl kplvec(minresidi) kplvec(minsosi) kplvec(minnsosi) kplvec(minnsosi2)];
kploi=[0  minresidi minsosi minnsosi minnsosi2];

parmso=[parms; kplcsfits(minresidi,:);  kplcsfits(minsosi,:);  kplcsfits(minnsosi,:); kplcsfits(minnsosi2,:)];

resido=[resid residvec(minresidi) residvec(minsosi) residvec(minnsosi) residvec(minnsosi2)];

if verbose
    f998h=figure(998);clf
    tiledlayout(4,1,'TileSpacing','Compact','Padding','Compact');
    nexttile
    plot(kplvec,residvec,'-',kpl,resid,'rx',kplvec(minresidi),minresid,'ro')
    %xlabel('k_{PL}')
    ylabel('Residual')
    title(sprintf('k_{PL} = %5.3f, %5.3f, %5.3f, %5.3f, %5.3f',kpl,kplvec(minresidi),kplvec(minsosi),...
        kplvec(minnsosi),kplvec(minnsosi2)));
    % nexttile % geometric mean of nuisance parameters
    % gmnp=geomean([kves;kecps;R1Lxs]); % doesn't work when any parm = 0;
    % %gmnp
    % semilogy(kplvec,gmnp)
    % ylabel('geomean(kve,kecp,R1Lx)')
    nexttile % sos of nuisance parms
    %sosnp=sqrt(kves.*kves+kecps.*kecps+R1Lxs.*R1Lxs);
    % size(kves)
    % sise(kecps)
    % size(R1Lxs)
    % size(sosnp)
    %[minsos minsosi]=min(sosnp);
    plot(kplvec,sosnp,'-',kplvec(minsosi),sosnp(minsosi),'rx');
    grid on
    ylabel('SOS')
    nexttile
    plot(kplvec,nsosnp,'-',kplvec(minnsosi),nsosnp(minnsosi),'rx');
    grid on
    ylabel('nSOS')
    nexttile
    plot(kplvec,nsosnp2,'-',kplvec(minnsosi2),nsosnp2(minnsosi2),'rx');
    grid on
    ylabel('nSOS2')
    xlabel('k_{PL}')    
    
    f999h=figure(999);clf
    tiledlayout(5,1,'TileSpacing','Compact','Padding','Compact');
   % nexttile
   % plot(kplvec,residvec,'-',kpl,resid,'rx',kplvec(minresidi),minresid,'ro')
   % ylabel('Residual')
    nexttile
    plot(kplvec,vbs,'-')
    title('Nuisance parameter trends');
    grid on
    ylabel('v_b (v/v)')
    nexttile
    plot(kplvec,kves,'-')
    grid on
    ylabel('K_{ve} (s^{-1})')
    nexttile
    plot(kplvec,vefs,'-')
    grid on
    ylabel('vef (v/v)')
    nexttile
    plot(kplvec,kecps,'-')
    grid on
    ylabel('K_{ecP} (s^{-1})')
    nexttile
    plot(kplvec,R1Lxs,'-')
    grid on
    ylabel('R_{1,Lx} (s^{-1})')
    xlabel('k_{PL} (s^{-1})')
    shg

end
