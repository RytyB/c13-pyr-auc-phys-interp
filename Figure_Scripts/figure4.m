%% Load results, if available
clear;
clc;
close all;

addpath('../Simulation_Functions');

%{
Figure 4

The point: Uncertainty in fitted parameters vary depending on whether
pyruvate spins are primarily vascular, extravascular/extracellular, or
intracellular.
%}

if exist('./ReproVsSNR.mat','file')
% if 0
    load('./ReproVsSNR.mat');
    fprintf('Found saved calculations - skip to last section to plot!\n')
else
    %% Set up model system:
    
    % Physics parameters
    base_knowns    = {'T1Pyr','T1Lac','T2spv','T2spe','T2spc','T2slc','Pe0','Pi0','Li0','FAp','FAl'};
    base_knownvals = [  43,     33,     0.1,   0.055   0.055,  0.033,   0,    0,    0,   20,   30];

    
    %Unknowns:
    fdve.fitvars={'kpl' 'VIFScale'};
    
    %Nuisance parameters; known or estimated eslewhere.
    % Note vef the extracellular fraction fraction:
    %      ve=(1-vb)*vef
    %      vc=1-vb-ve
    %physiological model parameters:
    fdve.knowns=   {'kve', 'kecp','vb',  'vef', 'R1Lx', base_knowns{:}} ;
    fdve.knownvals=[.0566, 0.0293, 0.1,   0.5,   0.0393, base_knownvals];

    %Describe acquisition scheme
    fdve.ntp=64; % Number of timepoints
    fdve.NSeg=1; % Segments per timepoint
    fdve.NFlips=(fdve.ntp)*(fdve.NSeg); % Total number of excitations
    
    %Describe temporal sampling scheme
    fdve.TE = 0.0219; % bbEPI; set TE=0 to ignore T2* differences
    TR=3; % seconds; constant repetition time
    fdve.TR=ones(1,fdve.NFlips)*TR;
    fdve.taxis=cumsum(fdve.TR)-fdve.TR(1);
    
    %Describe excitation scheme
    fdve.FlipAngle=[20 0;0 30]*ones(2,fdve.NFlips);
    
    %Describe vascular inmput function
    fdve.UseVIF=1;
    fdve.VIFP=gampdf(fdve.taxis-6,2.8,4.5);
    %fdv.VIFL=zeros(1,fdv.NFlips);
    
    %Placeholder for data to be fit
    fdve.data=ones(2,fdve.ntp);
    fdve.Name='Test Error in Parms, P2L3s';
    
    %Other miscellaneous
    fdve.verbose=0;
    
    %variables that are needed for slice profile
    fdve.SPBins=1;
    fdve.SPWeights=1;
    
    % Generate synthetic data:
    
    parms=[0.5371 1000]; %values for variables in fdv.fitvars
    parmse=parms;
    [EV, IV, vols]=P2L3s(parmse,fdve);
    
    fdve.data = [IV*vols(1) + EV(1,:)*vols(2) + EV(2,:)*vols(3);
                EV(3,:)*vols(3)];
    
    figure(61);clf
    set(gcf,'Name','EES Dominant')
    plot(fdve.taxis,fdve.data(1,:),'g-',...
        fdve.taxis,fdve.data(2,:),'b-');
    legend('Total Pyr','Total Lac')
    
    % confirm that parms lead to dominant signal from EES:
    [~, pvfse, ~, ~ ] = P2L3sAUC(parmse,fdve,fdve);
    fprintf('\nfdve: Pyr IV = %5.1f%%; EE = %5.1f%%; IC = %5.1f%%\n',pvfse*100)
    
    %% parms that reflect signal primarily from vasculature:
    
    fdvv=fdve;
    fdvv.knowns=   {'kve', 'kecp','vb',  'vef', 'R1Lx', base_knowns{:}};
    fdvv.knownvals=[.0113, 0.0253, 0.15,  0.5,  0.095,  base_knownvals];
    parmsv = [0.5645, 1000];
    
    [~, pvfsv, ~, ~ ] = P2L3sAUC(parmsv,fdvv,fdvv);
    fprintf('fdvv: Pyr IV = %5.1f%%; EE = %5.1f%%; IC = %5.1f%%\n',pvfsv*100)
    
    [EV, IV, vols]=P2L3s(parmsv,fdvv);
    
    fdvv.data = [IV*vols(1) + EV(1,:)*vols(2) + EV(2,:)*vols(3);
                EV(3,:)*vols(3)];
    
    figure(62);clf
    set(gcf,'Name','IV Dominant')
    plot(fdvv.taxis,fdvv.data(1,:),'g-',...
        fdvv.taxis,fdvv.data(2,:),'b-');
    legend('Total Pyr','Total Lac')
    
    
    %% parms that reflect signal primarily from inside cells:
    
    fdvc=fdve;
    fdvc.knowns=   {'kve', 'kecp','vb',  'vef', 'R1Lx', base_knowns{:}} ;
    fdvc.knownvals=[.1656, 0.121, 0.02,   0.2,   0.078, base_knownvals];
    parmsc = [0.2061, 1000]; % <-- lower kpl
    
    [~, pvfsc, ~, ~ ] = P2L3sAUC(parmsc,fdvc,fdvc);
    fprintf('fdvc: Pyr IV = %5.1f%%; EE = %5.1f%%; IC = %5.1f%%\n\n',pvfsc*100)
    
    [EV, IV, vols]=P2L3s(parmsc,fdvc);
    
    fdvc.data = [IV*vols(1) + EV(1,:)*vols(2) + EV(2,:)*vols(3);
                EV(3,:)*vols(3)];
    
    figure(63);clf
    set(gcf,'Name','IC Dominant')
    plot(fdvc.taxis,fdvc.data(1,:),'g-',...
        fdvc.taxis,fdvc.data(2,:),'b-');
    legend('Total Pyr','Total Lac')
    

    %% Set sweep parms and sweep!
    
    maxpyrsnr=10:10:100;    
    repspersnr=100;
    
    n_fits = 500;
    
    debug=0;
    jbopts=optimset('display','off');
    
    sigi=[4 45];
    
    kpls=zeros(3,length(maxpyrsnr),repspersnr);
    kves=zeros(3,length(maxpyrsnr),repspersnr);
    kecps=zeros(3,length(maxpyrsnr),repspersnr);
    ref_auc = zeros(3,1);
    aucs=zeros(3,length(maxpyrsnr),repspersnr);
    amps=zeros(3,length(maxpyrsnr),repspersnr);
    resids=zeros(3,length(maxpyrsnr),repspersnr);
    
    varvals=zeros(3,3);
    
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
        % add kve, kecp to fit parameters:
        fdvf.fitvars={fdvf.fitvars{:} 'kve' 'kecp'};
        parms=[parms fdvf.knownvals(1) fdvf.knownvals(2)];
        LB = [0 0 0 0];
        UB = [1 Inf 1 1];
        fdvf.knowns=fdvf.knowns(3:end);
        fdvf.knownvals=fdvf.knownvals(3:end);
        varvals(pvfi,:)=parms([3 4 1]);
        
        % ref_auc(pvfi) = Add this to make a reference line for fig4_pt5 

        fprintf('Testing %s, SNR = ',name);
        for ii=1:length(maxpyrsnr)
            fprintf('%d... ',maxpyrsnr(ii));
            for jj=1:repspersnr
                %add fresh noise each rep:
                fdvf.data=fdvk.data+(max(fdvk.data(1,:))/maxpyrsnr(ii))*randn(2,fdvf.ntp);
                %aucs(pvfi,ii,jj)=trapz(fdvf.data(2,sigi(1):sigi(2)).*fdvf.TR(sigi(1):sigi(2)))./trapz(fdvf.data(1,sigi(1):sigi(2)).*fdvf.TR(sigi(1):sigi(2)));
                %run fits:
                %bestresid=Inf;
                %bestfits = zeros(1,length(parms));
                fitparms=zeros(n_fits,length(parms));
                fitresids=zeros(1,n_fits);
                parfor kk=1:n_fits
                    %nrkpl=zeros(1,n_fits);
                    %ii
                    Guess=([1 5000 1 1]).*rand(1,length(parms));
                    [tmpparms,tmpresid] = lsqnonlin(@(x) P2L3sErr(x,fdvf),Guess,LB,UB,jbopts);
                    fitparms(kk,:)=tmpparms;
                    fitresids(kk)=tmpresid;
                    %if resid<bestresid
                    %    bestresid=resid;
                    %    bestfits=fits;
                    %end
                end
                [bestresid, bfi]=min(fitresids);
                aucs(pvfi,ii,jj) = P2L3sAUC(fitparms(bfi,:),fdvf,fdvf);
                kpls(pvfi,ii,jj)=fitparms(bfi,1);%bestfits(1);
                amps(pvfi,ii,jj)=fitparms(bfi,2);%bestfits(2);
                kves(pvfi,ii,jj)=fitparms(bfi,3);%bestfits(3);
                kecps(pvfi,ii,jj)=fitparms(bfi,4);%bestfits(4);
                resids(pvfi,ii,jj)=bestresid;
                if debug
                    [EV, IV, vols]=P2L3s(bestfits,fdvf);
                    fdata = [IV*vols(1) + EV(1,:)*vols(2) + EV(2,:)*vols(3);
                        EV(3,:)*vols(3)];
                    figure(99);clf
                    plot(fdvk.taxis,fdvf.data(1,:),'gx',...
                        fdvk.taxis,fdvf.data(2,:),'bx',...
                        fdvk.taxis,fdata(1,:),'g-',...
                        fdvk.taxis,fdata(2,:),'b-');
                    legend('Total Pyr','Total Lac')
                    title(sprintf('k_{PL} = %5.3f',bestfits(1)));
                    xlabel('Time (s)')
                    pause
                end
                %fprintf('Kpla in = %5.3f; Kpla out = %5.3f\n',parms(1),bestfits(1));
            end % snr loop
        end % repspersnr loop
        fprintf('Done!\n')
    end % pvfi loop
    
    if ~exist('../Data/ReproVsSNR.mat','file')
        save('../Data/ReproVsSNR.mat')
    end

end  % End if data exists


%% Now plot!
fig4 = figure('Name', 'Figure 4');
tiledlayout(3,3,'TileSpacing','Compact','Padding','Compact');
fontsize = 12;
linewidth = 1.25;

fig4.Units = 'inches';
fig4.Position = [7.0729, 1.9479, 6.9, 6.9];

for pvfi=1:3
    %Plot curves for  kve:
    nexttile(pvfi)
    % yyaxis left
    errorbar(maxpyrsnr,mean(squeeze(kves(pvfi,:,:))'),std(squeeze(kves(pvfi,:,:))'),...
        'LineWidth',linewidth);
    refline(0,varvals(pvfi,1));
    set(gca, 'FontSize',fontsize*0.8);
    if pvfi==1      
        ylabel('k_{ve} (s^{-1})','FontSize',fontsize, 'FontWeight','bold');
        title('Pyr_{iv} Dominant','FontSize',fontsize);
    elseif pvfi==2
        title('Pyr_{ee} Dominant','FontSize',fontsize);
    elseif pvfi==3
        title('Pyr_c Dominant','FontSize',fontsize);
    end
    axis tight
    grid on

        
    % Plot curves for kec:
    nexttile(3+pvfi)
    % yyaxis left
    errorbar(maxpyrsnr,mean(squeeze(kecps(pvfi,:,:))'),std(squeeze(kecps(pvfi,:,:))'),...
        'LineWidth',linewidth);
    refline(0,varvals(pvfi,2));
    set(gca, 'FontSize',fontsize*0.8);
    if pvfi==1
        ylabel('k_{ecP} (s^{-1})','FontSize',fontsize, 'FontWeight','bold')
    end
    axis tight
    % yyaxis right
    % plot(maxpyrsnr,std(squeeze(kecps(pvfi,:,:))')./mean(squeeze(kecps(pvfi,:,:))'),'linewidth',2)
    % axis tight
    % if pvfi==3
    %     ylabel('Coefficient of Variation')
    % end
    grid on


    %plot curves for kpl:
    nexttile(6+pvfi)
    % yyaxis left
    errorbar(maxpyrsnr,mean(squeeze(kpls(pvfi,:,:))'),std(squeeze(kpls(pvfi,:,:))'),...
        'LineWidth',linewidth);
    refline(0,varvals(pvfi,3));
    set(gca, 'FontSize',fontsize*0.8);
    if pvfi==1
        ylabel('k_{PL} (s^{-1})','FontSize',fontsize, 'FontWeight','bold');
    end
    axis tight;
    % yyaxis right
    % plot(maxpyrsnr,std(squeeze(kpls(pvfi,:,:))')./mean(squeeze(kpls(pvfi,:,:))'),'linewidth',2)
    % axis tight
    % if pvfi==3
    %     ylabel('Coefficient of Variation')
    % end
    grid on;
    xlabel('Max Pyr SNR','FontSize',fontsize);
end % pvfi loop

% Change line colors to match other figures
red = [0.735, 0.078, 0.184];
green = [0.194, 0.529, 0.025];
yellow = [0.929, 0.694, 0.125]; 

for ii=1:3:9
    nexttile(ii);
    iv_plot = gca;
    iv_plot.Children(1).Color = 'k';
    iv_plot.Children(2).Color = red;
    
    nexttile(ii+1);
    ee_plot = gca;
    ee_plot.Children(1).Color = 'k';
    ee_plot.Children(2).Color = yellow;

    nexttile(ii+2);
    c_plot = gca;
    c_plot.Children(1).Color = 'k'; 
    c_plot.Children(2).Color = green; 
end


%% Save pretty picture and data from simulation
% Export as MATLAB fig file for later editing
savefig(fig4, '../Media/figure4.fig');

% Export as vector image for publication
exportgraphics(fig4, '../Media/figure4.eps',...
    'Resolution', 600,...
    'BackgroundColor', 'white',...
    'ContentType', 'vector');

% Export as pixel image for preview
exportgraphics(fig4, '../Media/figure4.png',...
    'Resolution', 600,...
    'BackgroundColor', 'white');


%% Figure 4.5 ??
fig4_pt5 = figure('Name','AUCr??');
fig4_pt5.Units = 'inches';
fig4_pt5.Position = [0.0833, 6.1667, 6.9, 2.6771];

tiledlayout(1,3, 'TileSpacing','compact');

nexttile;
errorbar(maxpyrsnr,mean(squeeze(aucs(1,:,:))'),std(squeeze(aucs(1,:,:))'),...
        'Color',red, 'LineWidth',linewidth);
set(gca, 'FontSize',fontsize*0.8);
xlabel('Max Pyr SNR', 'FontSize',fontsize);
ylabel('Calcd. AUC Ratio', 'FontSize',fontsize);
title('Pyr_{iv} Dominant', 'FontSize',fontsize);
grid on;
axis tight;
xlim([0,maxpyrsnr(end)]); 

nexttile;
errorbar(maxpyrsnr,mean(squeeze(aucs(2,:,:))'),std(squeeze(aucs(2,:,:))'),...
        'Color',yellow, 'LineWidth',linewidth);
set(gca, 'FontSize',fontsize*0.8);
xlabel('Max Pyr SNR', 'FontSize',fontsize);
title('Pyr_{ee} Dominant', 'FontSize',fontsize);
grid on;
axis tight;
xlim([0,maxpyrsnr(end)]);

nexttile;
errorbar(maxpyrsnr,mean(squeeze(aucs(3,:,:))'),std(squeeze(aucs(3,:,:))'),...
        'Color',green, 'LineWidth',linewidth);
set(gca, 'FontSize',fontsize*0.8);
xlabel('Max Pyr SNR', 'FontSize',fontsize);
title('Pyr_{c} Dominant', 'FontSize',fontsize);
grid on;
axis tight;
xlim([0,maxpyrsnr(end)]);