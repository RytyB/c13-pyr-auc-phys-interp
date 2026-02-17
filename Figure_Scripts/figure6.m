%% Prepare Workspace:
clear
close all

infn='PVFsVsSNRParmSweep.20250909.mat';

addpath('../Simulation_Functions');

%{
Figure 6

The point: Even if you can't get reproducible parameter fits, 
you can do pretty well at figuring out which of the compartments 
the pyruvate comes from and thus figure out the rate limiting step.  
%}


%% Create model systems:

% Unknowns:
fdve.fitvars={'kpl' 'VIFScale' 'kve'  'kecp' 'vb' 'vef' 'R1Lx'}; % Parameters for fitting
parmse =     [0.5371  1000     0.0566 0.0293  0.1  0.5  0.0393];
% Nuisance parameters; known or estimated eslewhere.
% Note vef the extracellular fraction fraction:
%      ve=(1-vb)*vef
%      vc=1-vb-ve

% Physiological model parameters:
% fdve.knowns=  {'R1Lx'} ;
% fdve.knownvals=[0.0];

% Physics parameters:
fdve.knowns=   {'T1Pyr','T1Lac','T2spv','T2spe','T2spc','T2slc','Pe0','Pi0','Li0','FAp','FAl'};
fdve.knownvals=[43,     33,     0.1,    0.055   0.055,  0.033,  0,    0,    0,    20,   30];

% Define vascular input function (VIF)
fdve.knowns=   {fdve.knowns{:},'Gam1','Gam2','tdel',};
fdve.knownvals=[fdve.knownvals, 2.8,   4.5,   6];

% Describe acquisition scheme
fdve.ntp=64; % Number of timepoints
fdve.NSeg=1; % Segments per timepoint
fdve.NFlips=(fdve.ntp)*(fdve.NSeg); % Total number of excitations

% Describe temporal sampling scheme
fdve.TE = 0.0; %0.0219; % Echo time, bbEPI; set TE=0 to ignore T2* differences
TR=2; % Repitition time (seconds) & constant
fdve.TR=ones(1,fdve.NFlips)*TR;
fdve.taxis=cumsum(fdve.TR)-fdve.TR(1);

% Describe excitation scheme
fdve.FlipAngle=[20 0;0 30]*ones(2,fdve.NFlips); % Flip angle matrix

% Describe vascular input function usage
fdve.UseVIF=1; % Enable VIF usage
fdve.VIFP=gampdf(fdve.taxis-6,2.8,4.5); % Define based on gamma distribution

% Placeholder for synthetic data
fdve.data=ones(2,fdve.ntp); % Initialize data matrix
fdve.Name='Test Error in Parms, P2L3s'; % Experiment name

% Miscellaneous parameters
fdve.verbose=0; % Disable verbose output
fdve.SPBins=1; % Assume ideal slice profile
fdve.SPWeights=1;

% Generate synthetic data:

[EV, IV, vols]=P2L3s(parmse,fdve); % Run model simulation
fdve.data = [IV*vols(1) + EV(1,:)*vols(2) + EV(2,:)*vols(3);
            EV(3,:)*vols(3)];
% Plot generated data
figure(71);clf
set(gcf,'Name','EES Dominant')
plot(fdve.taxis,fdve.data(1,:),'g-',...
    fdve.taxis,fdve.data(2,:),'b-');
legend('Total Pyr','Total Lac')
xlabel('Time (s)');
ylabel('Amplitude (arb)')

% These parms lead to dominant signal from EES:
 [~, ~, pvfse,~] = P2L3sAUC(parmse,fdve);
fprintf('\nfdve: Pyr IV = %5.1f%%; EE = %5.1f%%; IC = %5.1f%%\n',pvfse*100)

%% Parms that reflect signal primarily from vasculature:
fdvv=fdve;
%      {'kpl' 'VIFScale' 'kve'  'kecp' 'vb' 'vef' 'R1Lx'}
parmsv=[0.5645  1000    0.0113  0.0253  0.15  0.5  0.0095];
 [~, ~, pvfsv, ~ ] = P2L3sAUC(parmsv,fdvv);
fprintf('fdvv: Pyr IV = %5.1f%%; EE = %5.1f%%; IC = %5.1f%%\n',pvfsv*100)
[EV, IV, vols]=P2L3s(parmsv,fdvv);
fdvv.data = [IV*vols(1) + EV(1,:)*vols(2) + EV(2,:)*vols(3);
            EV(3,:)*vols(3)];
% Plot generated data 
figure(72);clf
set(gcf,'Name','IV Dominant')
plot(fdvv.taxis,fdvv.data(1,:),'g-',...
    fdvv.taxis,fdvv.data(2,:),'b-');
legend('Total Pyr','Total Lac')
xlabel('Time (s)');
ylabel('Amplitude (arb)')

%% Parms reflect signal primarily from inside cells:
fdvc=fdve;
%      {'kpl' 'VIFScale' 'kve'  'kecp' 'vb' 'vef' 'R1Lx'}
parmsc=[0.2061  1000    0.1656  0.1211  0.02 0.2  0.078];
%parmsc=parms;
parmsc(1)= 0.2; % lower kpl
 [~, ~, pvfsc,~] = P2L3sAUC(parmsc,fdvc);
fprintf('fdvc: Pyr IV = %5.1f%%; EE = %5.1f%%; IC = %5.1f%%\n\n',pvfsc*100)
[EV, IV, vols]=P2L3s(parmsc,fdvc);
fdvc.data = [IV*vols(1) + EV(1,:)*vols(2) + EV(2,:)*vols(3);
            EV(3,:)*vols(3)];
% Plot generated data
figure(73);clf
set(gcf,'Name','IC Dominant')
plot(fdvc.taxis,fdvc.data(1,:),'g-',...
    fdvc.taxis,fdvc.data(2,:),'b-');
legend('Total Pyr','Total Lac')
xlabel('Time (s)');
ylabel('Amplitude (arb)')

%% Range of max pyr SNR values that will be tested:
    
if exist(infn,'file')
    fprintf('Found saved data! Loading in data...\n')
    load(infn)
else % run simulations!
    maxpyrsnr=10:10:100; % Range of SNR values
    repspersnr=50; % Number of repitions per SNR
    
    kplveclen=30;
    kplvec=logspace(log10(0.01),log10(1.5),kplveclen);
    
    % nruns=2;
    nruns=5000; % Number of runs per rep
    
    debug=0;
    jbopts=optimset('display','off');
    sigi=[4 45];
    
    % Initialize variables for storage
    % test 5 wasy of estimating kpl:
    kpls=zeros(3,5,length(maxpyrsnr),repspersnr);
    kves=kpls;
    aucs=kpls;
    amps=kpls;
    resids=kpls;
    kecps=kpls;
    pvfs=zeros(3,5,length(maxpyrsnr),repspersnr,3);
    %varvals=zeros(3,5,3);
    
    for pvfi=1:3
        switch pvfi
            case 1 %pyr signal mostly from vasculature
                fdvk=fdvv;
                parms=parmsv;
                name='PyrV';
            case 2 %pyr signal mostly from EES
                fdvk=fdve;
                parms=parmse;
                name='PyrE';
            case 3 %pyr signal mostly from intracellular space
                fdvk=fdvc;
                parms=parmsc;
                name='PyrC';
        end
    
        fdvf=fdvk;
        fdvf.Name=name;
        % move kpl from unknowns for the sweep:
        fdvf.knowns=   {fdvf.knowns{:},'kpl'};
        fdvf.knownvals=[fdvf.knownvals parms(1)];
        fdvf.fitvars=fdvf.fitvars(2:end);
    
        
        % Upper and lower bounds
        %LB   =[0 0   0 0 0 (eps) 0];
        %UB   =[1 Inf 1 1 1 (1-eps) 1];
        LB   =[ 0          0      0      0    (eps)    0];
        UB   =[ Inf        0.2    0.2      1    (1-eps)  0.1];
        LGB  =[ 0          1e-4   1e-4   0    (eps)    0];
        UGB  =[ 0          0.2    0.2    1    (1-eps)  0.1];
        %varvals(pvfi,:)=parms([3 4 1]);
        varvals=[pvfsv ; pvfse ; pvfsc];
    
        fprintf('Testing %s, SNR = ',name);
        for ii=1:length(maxpyrsnr)
            fprintf('%d... ',maxpyrsnr(ii));
            for jj=1:repspersnr
                fdvf.Name=sprintf('%s - %d/%d, %d/%d',name,ii,length(maxpyrsnr),jj,repspersnr);
                %add fresh noise each rep:
                fdvf.data=fdvk.data+(max(fdvk.data(1,:))/maxpyrsnr(ii))*randn(2,fdvf.ntp);
                aucs(ii,jj)=trapz(fdvf.data(2,sigi(1):sigi(2)).*fdvf.TR(sigi(1):sigi(2)))./trapz(fdvf.data(1,sigi(1):sigi(2)).*fdvf.TR(sigi(1):sigi(2)));
                
                [okpl, oparms, resid, residvec]=P2L3sKPLSweep(fdvf,kplvec,nruns,LB,UB,LGB,UGB);
                for kpli = 1:5
                    %[~,~,pvf,~]=P2L3sAUC(bestfits,fdvf);%best fit for loop of runs
                    fdvf.knownvals(end)=okpl(kpli); % sub in kpl from sweep
                    [~,~,pvf,~]=P2L3sAUC(oparms(kpli,:),fdvf);
                    pvfs(pvfi,kpli,ii,jj,:)=pvf;
                end
                
            end % repssnr loop
        end % snr loop
        fprintf('Done!\n')
    end % pvfi loop
    
    tmp=datetime('today');
    ofn=sprintf('%s.%d%02d%02d.mat',mfilename,tmp.Year,tmp.Month,tmp.Day);
    copyfile([mfilename '.m'],sprintf('%s.%d%02d%02d.m',mfilename,tmp.Year,tmp.Month,tmp.Day));

    if ~exist(ofn,'file')
        save(ofn)
    end
end
%% Now plot!

kpld={'1pct', 'BestFit', 'SOS', 'nSOS', 'nSOS2'};


%% what fraction of these analyses identified the correct max pool:

tp=zeros(3,5,length(maxpyrsnr)); % true positives
tn=tp;                           % true negatives
fn=tp;                           % false negatives
fp=tp;                           % false positives
pc=zeros(3,5,3,length(maxpyrsnr)); %positive calls

% Sensitivity : TP/(TP+FN), noting TP+FN = repspersnr
% Specificity : TN/(TN+FP), noting TN+FP = 2*repspersnr

for pvfi=1:3
    for kpldi=1:5
        for ii=1:length(maxpyrsnr)
            for jj=1:repspersnr
                [~,maxi]=max(pvfs(pvfi,kpldi,ii,jj,:));
                %maxi
                if maxi==pvfi
                    tp(pvfi,kpldi,ii)=tp(pvfi,kpldi,ii)+1;
                    if pvfi~=1
                        tn(1,kpldi,ii)=tn(1,kpldi,ii)+1;
                    end
                    if pvfi~=2
                        tn(2,kpldi,ii)=tn(2,kpldi,ii)+1;
                    end
                    if pvfi~=3
                        tn(3,kpldi,ii)=tn(3,kpldi,ii)+1;
                    end
                else
                    fn(pvfi,kpldi,ii)=fn(pvfi,kpldi,ii)+1;
                    fp(maxi,kpldi,ii)=fp(maxi,kpldi,ii)+1;
                end
                pc(pvfi,kpldi,maxi,ii)=pc(pvfi,kpldi,maxi,ii)+1;
            end
        end
    end
end


f177h=figure(177);
fontsize = 14;
f177h.Units = 'inches';
f177h.Position=[2.2083, 1.1417, 4.2333, 6.9500];
tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact');
for kpldi=5
    for pvfi=1:3
        %sens=tp(pvfi,kpldi,:)./(tp(pvfi,kpldi,:)+fn(pvfi,kpldi,:));
        sens=squeeze(tp(pvfi,kpldi,:)./(tp(pvfi,kpldi,:)+fn(pvfi,kpldi,:)));
        spec=squeeze(tn(pvfi,kpldi,:)./(tn(pvfi,kpldi,:)+fp(pvfi,kpldi,:)));
        nexttile
        %yyaxis left
        tmph=plot(maxpyrsnr,100*squeeze(sens),maxpyrsnr,100*squeeze(spec),'LineWidth',2);
        set(gca,'FontSize',11);
        axis([0 max(maxpyrsnr) 0 100]);
        grid on        
        switch pvfi
            case 1
                title('Pyr_{iv} Dominant','FontSize',fontsize);
                tmph=legend('Sensitivity','Specificity');
                tmph.FontSize=12;
                fprintf('PyrIV Avg Sensitivity = %5.3f, Avg Specificity = %5.3f\n',mean(sens),mean(spec))
            case 2
                title('Pyr_{ee} Dominant','FontSize',fontsize);
                fprintf('PyrEE Avg Sensitivity = %5.3f, Avg Specificity = %5.3f\n',mean(sens),mean(spec))
            case 3
                title('Pyr_c Dominant','FontSize',fontsize);
                fprintf('PyrIC Avg Sensitivity = %5.3f, Avg Specificity = %5.3f\n',mean(sens),mean(spec))
        end
        if pvfi==2
            ylabel('Percentage','FontSize',fontsize);
        end
        % yyaxis right
        % plot(maxpyrsnr,100*squeeze(spec),'LineWidth',2);
        % axis([0 max(maxpyrsnr) 0 100]);
        % if pvfi==3
        %     ylabel('Specificity');
        % end

        ylim([40,100]);
    end
    xlabel('Max Pyr SNR','FontSize',fontsize);

end
ofigfn=sprintf('%s-FCorrect-SOS2.%d%02d%02d.png',mfilename,tmp.Year,tmp.Month,tmp.Day);
if ~exist(ofigfn,'file')&&(~strncmp(ofigfn,'Live',4))
    exportgraphics(gcf,ofigfn)
end

