function [Mtev, Mtiv, vols, Mzev, Mziv] = P2L3s(vars,fdv)
%Format: function [Mxyev, Mxyiv, vols, Mev, Miv] = P2L3sSP(vars,fdv)
%
%input:
    % vars:                 values for parms have been or will be fit
    % fdv = "forward variables" containing info about the measurement
    %    fdv.fitvars        Names of parms that need to be fit
    %    fdv.knowns         Names of parms that are known
    %    fdv.knownvals      Values of known parameters
    %    fdv.FlipAngles     Excitation angle list [2 x NFlips] in degrees
    %    fdv.TR:            Repetition time [2 x NFlips]
    %    fdv.taxis:         time axis for the data [1 x NFlips]
    %    fdv.data           [2 x NTP] observed [pyr;lac] data
    %    fdv.UseVIF         =1 if VIFs have been measured
    %    fdv.VIFP           **Mz** VIF for pyruvate [1 x NFlips]
    %    fdv.VIFL           **Mz** VIF for lactate [1 x NFlips] (usually zeros)
    %    fdv.verbose        =1 to plot results of this script
%    
%output:
    % Mxyev                 [3 x NTP] extravascular observed [PyrEE;PyrIC;LacIC] 
    % Mxyiv                 [1 x NTP] IV observed [PyrIV]
    % vols                  vascular blood vol fractions [vb ve vc]            
    % Mzev                  [3 x NTP] extravascular longitudinal [PyrEE;PyrIC;LacIC]
    % Mziv                  [1 x NTP] IV longitudinal [PyrIV];

%
% Note Mz reflects magnetization before excitation pulse 
%  while Mxy reflects transverse magnetization due to pulse.
%
% Note total observed signal is 
% [Pyr;Lac] = vb*Mxyiv + ve*Mxyev(1:2) + vc*Mxyev(3:4)
%           = vb*Mxyiv + (1-vb)*vef*Mxyev(1:2) + (1-(1-vb)*vef)*Mxyev(3:4)

%    eps=1e-6;
       
%Unpack variables:
for ii=1:length(vars)
    x.(fdv.fitvars{ii})=vars(ii);
end
for ii=1:length(fdv.knowns)
    x.(fdv.knowns{ii})=fdv.knownvals(ii);
end
if isfield(x,{'FA'}) % sub fit values for fdv values
    fdv.FlipAngle=x.FA*ones(size(fdv.FlipAngle));
end
if isfield(x,{'FAp'}) % sub fit values for fdv values
    fdv.FlipAngle=[x.FAp 0 ; 0 x.FAl]*ones(size(fdv.FlipAngle));
end

% vb=x.vb;
% ve=(1-vb)*x.vef;
% vc=1-vb-ve;
% vols=[vb ve vc];
vols=[x.vb (1-x.vb)*x.vef (1 - x.vb - x.vef + x.vb*x.vef)];
R1Pyr=1/x.T1Pyr;
R1Lac=1/x.T1Lac;
kvedve=x.kve/vols(2);
kecpdve=x.kecp/vols(2);
kecpdvc=x.kecp/vols(3);

%Assume for now that the VIF is given:
if fdv.UseVIF
    Mziv=x.VIFScale*fdv.VIFP;
else 
    Mziv=x.VIFScale*gampdf(fdv.taxis-x.tdel,x.Gam1,x.Gam2);
end

% %IV magnetization at each observation:
% Mtiv=Mziv.*sind(fdv.FlipAngle(1,:))*exp(-fdv.TE/x.T2spv);

%Initialize Variables:
Pze=zeros(1,fdv.ntp); % Longitudinal pyruvate in extravascular/extracellular space 
Pzc=zeros(1,fdv.ntp); % Longitudinal pyr in intracellular space
Lzc=zeros(1,fdv.ntp); % Longitudinal lactate in intracellular space
Ptv=zeros(1,fdv.ntp); % Transverse pyr in vascular space
Pte=zeros(1,fdv.ntp); % Transverse pyr in EES
Ptc=zeros(1,fdv.ntp); % Transverse pyr, intracellular
Ltc=zeros(1,fdv.ntp); % Transverse lac, intracellular
%Initialize return variables:
Mtiv = zeros(1,fdv.ntp); 
Mzev = [Pze; Pzc; Lzc]; 
Mtev = [Pte; Ptc; Ltc]; 
%
% Here the equations are decoupled so that for each phys/chem
% compartment we have:
%
% Diff Eq: y'(t) = Ay(t) + ff(t)
% Sol'n:   y(t) = y(0)*exp(A*t) + exp(A*t)*integral(0,t: exp(-A*tau)*ff(tau) dtau)
%

%attenuation coefficient for extravascular/extracellular pyr:
ape = -(kvedve+kecpdve+R1Pyr);
%attenuation coefficient for intracellular pyr:
%apc = -(kecpdvc+x.kpl+R1Pyr);
apc = -(x.kpl+R1Pyr);
%attenuation coefficient for intracellular lac, including extra loss term:
alc = -(R1Lac+x.R1Lx);


%calculate as much as possible outside of the for loop:

%VIF - assume linear connection between observations separated by TR: 
% Pzb(t)=mt+b        
TR=fdv.TR(1:end-1); %Note first TR is delay between first and second observations
vifm=(Mziv(2:end)-Mziv(1:end-1))./TR;
vifb=Mziv(1:end-1);
%Terms that only depend on the VIF:
PE1 = -kvedve*(vifb/ape+vifm/ape/ape);
PE2 = -kvedve*vifm/ape;
PC1 = -kecpdvc*(PE1/apc+PE2/apc/apc);
PC2 = -kecpdvc*PE2/apc;  
LC1 = -x.kpl*(PC1/alc+PC2/alc/alc);
LC2 = -x.kpl*PC2/alc;

%Loop through bins of slice profile:
for zz=1:length(fdv.SPBins)
    %initial conditions:
    Pe0=x.Pe0; 
    Pi0=x.Pi0; 
    Li0=x.Li0;
    % Mtiv=Mziv.*sind(fdv.FlipAngle(1,:))*exp(-fdv.TE/x.T2spv);
    %Calculate dynamic signal evolution:
    for ii=1:(fdv.ntp-1)
        Ptv(ii)=Mziv(ii)*sind(fdv.FlipAngle(1,ii)*fdv.SPBins(zz))*exp(-fdv.TE/x.T2spv);
        Pze(ii)=Pe0;
        Pte(ii)=Pe0*sind(fdv.FlipAngle(1,ii)*fdv.SPBins(zz))*exp(-fdv.TE/x.T2spe);
        Pzc(ii)=Pi0;
        Ptc(ii)=Pi0*sind(fdv.FlipAngle(1,ii)*fdv.SPBins(zz))*exp(-fdv.TE/x.T2spc);
        Lzc(ii)=Li0;
        Ltc(ii)=Li0*sind(fdv.FlipAngle(2,ii)*fdv.SPBins(zz))*exp(-fdv.TE/x.T2slc);
        %Calculate magnetization at the end of this TR = IC for next TR:
        % Extravascular/extracellular pyruvate:
        PE3 = Pe0*cosd(fdv.FlipAngle(1,ii)*fdv.SPBins(zz)) - PE1(ii);
        Pe0 = PE1(ii) + PE2(ii)*TR(ii) + PE3*exp(ape*TR(ii));
        % Intracellular Pyruvate:
        PC4 = kecpdvc*PE3/(ape-apc);
        PC3 = Pi0*cosd(fdv.FlipAngle(1,ii)*fdv.SPBins(zz)) - PC1(ii) - PC4;
        Pi0 = PC1(ii) + PC2(ii)*TR(ii) + PC3*exp(apc*TR(ii)) + PC4*exp(ape*TR(ii));
        % Intracellular Lactate:
        LC4 = x.kpl*PC3/(apc-alc);
        LC5 = x.kpl*PC4/(ape-alc);
        LC3 = Li0*cosd(fdv.FlipAngle(2,ii)*fdv.SPBins(zz)) - LC1(ii) - LC4 - LC5;
        Li0 = LC1(ii) + LC2(ii)*TR(ii) + LC3*exp(alc*TR(ii)) + LC4*exp(apc*TR(ii)) + LC5*exp(ape*TR(ii));
    end
    Ptv(end)=Mziv(end)*sind(fdv.FlipAngle(1,end)*fdv.SPBins(zz))*exp(-fdv.TE/x.T2spv);
    Pze(end)=Pe0;
    Pte(end)=Pe0*sind(fdv.FlipAngle(1,end)*fdv.SPBins(zz));
    Pzc(end)=Pi0;
    Ptc(end)=Pi0*sind(fdv.FlipAngle(1,end)*fdv.SPBins(zz));
    Lzc(end)=Li0;
    Ltc(end)=Li0*sind(fdv.FlipAngle(2,end)*fdv.SPBins(zz));
    
    %Fill return variables:
    Mtiv = Mtiv + Ptv*fdv.SPWeights(zz);
    Mzev = Mzev + [Pze; Pzc; Lzc]*fdv.SPWeights(zz); % Longitudinal magnetization
    Mtev = Mtev + [Pte; Ptc; Ltc]*fdv.SPWeights(zz); % Transverse magnetization
end % slice profile bin loop

if fdv.verbose
    figure(99)
    plot(fdv.taxis(1:end),Mtiv(1,:)*vols(1),'r:',...
        fdv.taxis(1:end),Mtev(1,:)*vols(2),'g:',...
        fdv.taxis(1:end),Mtev(2,:)*vols(3),'g-',...
        fdv.taxis(1:end),Mtev(3,:)*vols(3),'b-')
    legend('PyrIV','PyrEV','PyrIC','LacIC')
end    