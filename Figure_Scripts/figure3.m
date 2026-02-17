%% Set up model system:
clear; 
close all; 
clc;

addpath('../Simulation_Functions');

%{
Figure 3

The point: The AUC ratio changes more for changes in underlying kinetic
parameters depending on whether pyruvate is primarily vascular,
extravascular/extracellular, or intracellular.
%}


%Unknowns:
fdve.fitvars={'kpl' 'VIFScale'};

%physics parameters:
base_knowns=   {'T1Pyr','T1Lac','T2spv','T2spe','T2spc','T2slc','Pe0','Pi0','Li0','FAp','FAl'};
base_knownvals=[  43,     33,     0.1,    0.055   0.055,  0.033,  0,    0,    0,    20,   30];

%Nuisance parameters; known or estimated eslewhere.
% Note vef the extracellular fraction fraction:
%      ve=(1-vb)*vef
%      vc=1-vb-ve
%% parms that reflect signal primarily from extracellular/extravascular space
fdve.knowns=   {base_knowns{:}, 'kve', 'kecp','vb',  'vef', 'R1Lx'} ;
fdve.knownvals=[base_knownvals, .0566, 0.0293, 0.1,   0.5,   0.0393];
%parameters for the VIF (soon will be deprecated)
fdve.knowns=   {fdve.knowns{:},'Gam1','Gam2','tdel',};
fdve.knownvals=[fdve.knownvals, 2.8,   4.5,   6];


%fdv.knowns=  {'kve','kecp','vb', 'vef', 'T1Pyr','T1Lac','T2spv','R1Lx','Gam1','Gam2','tdel','Pe0','Pi0','Li0'};
%fdv.knownvals=[0.02, 0.02, 0.05,  0.5,   43,     33,     0.1,    0,    2.8,   4.5,   10,    0,     0,    0];

%Describe acquisition scheme
fdve.ntp=64; % Number of timepoints
fdve.NSeg=1; % Segments per timepoint
fdve.NFlips=(fdve.ntp)*(fdve.NSeg); % Total number of excitations

%Describe temporal sampling scheme
fdve.TE = 0.0219; % bbEPI; set TE=0 to ignore T2* differences
% fdve.TE = 0;
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
[EV, IV, vols]=P2L3s(parms,fdve);

%fdv.data=IV(1:2,:)*vols(1) + EV(1:2,:)*vols(2) + EV(3:4,:)*vols(3);
fdve.data = [IV*vols(1) + EV(1,:)*vols(2) + EV(2,:)*vols(3);
            EV(3,:)*vols(3)];
%resid=P2LBv4Err(parms,fdv);

% confirm that parms lead to dominant signal from EES:

[~, pvfse, ~, ~] = P2L3sAUC(parms,fdve,fdve);
fprintf('\nEE Dominant: Pyr IV = %5.1f%%; EE = %5.1f%%; IC = %5.1f%%\n',pvfse*100)

%% parms that reflect signal primarily from vasculature:

fdvv=fdve;
fdvv.knowns=   {base_knowns{:}, 'kve', 'kecp','vb',  'vef', 'R1Lx'};
fdvv.knownvals=[base_knownvals, .0113, 0.0253, 0.15,  0.5,   0.095];
parms = [0.5645, 1000];

[~, pvfsv, ~, ~ ] = P2L3sAUC(parms,fdvv,fdvv);
fprintf('IV Dominant: Pyr IV = %5.1f%%; EE = %5.1f%%; IC = %5.1f%%\n',pvfsv*100)

[EV, IV, vols]=P2L3s(parms,fdvv);

%fdv.data=IV(1:2,:)*vols(1) + EV(1:2,:)*vols(2) + EV(3:4,:)*vols(3);
fdvv.data = [IV*vols(1) + EV(1,:)*vols(2) + EV(2,:)*vols(3);
            EV(3,:)*vols(3)];
%resid=P2LBv4Err(parms,fdv);




%% parms that reflect signal primarily from inside cells:

fdvc=fdve;
fdvc.knowns=   {base_knowns{:}, 'kve', 'kecp','vb',  'vef', 'R1Lx'} ;
fdvc.knownvals=[base_knownvals, .1656, 0.121, 0.02,   0.2,   0.078];
parms = [0.2061, 1000]; % <-- lower kpl

[~, pvfsc, ~, ~ ] = P2L3sAUC(parms,fdvc,fdvc);
fprintf(' C Dominant: Pyr IV = %5.1f%%; EE = %5.1f%%; IC = %5.1f%%\n\n',pvfsc*100)

[EV, IV, vols]=P2L3s(parms,fdvc);

%fdv.data=IV(1:2,:)*vols(1) + EV(1:2,:)*vols(2) + EV(3:4,:)*vols(3);
fdvc.data = [IV*vols(1) + EV(1,:)*vols(2) + EV(2,:)*vols(3);
            EV(3,:)*vols(3)];
%resid=P2LBv4Err(parms,fdv);


%% Calculate changes in AUC vs changes in rate-limiting parms:

%testparms={'vb','kve','kecp','kpl'};
testparms = {'kve', 'kecp', 'kpl'}; 
nefpoints=21;
ef=linspace(-0.25, 0.25, nefpoints);
aucs=zeros(3,length(testparms),nefpoints);

for ii=1:3 
    switch ii % switch between parameter sets
        case 1
            fdvr=fdvv;
        case 2
            fdvr=fdve;
        case 3
            fdvr=fdvc;
    end
    for jj=1:length(testparms)
        fdvt=fdvr;
        parmst=parms;
        if strcmp(testparms{jj},'kpl')
            parmi = find(strcmp(fdvt.fitvars,'kpl'));
        else
            kpi = find(strcmp(fdvt.knowns,testparms{jj}));
        end
        for kk=1:nefpoints
            if strcmp(testparms{jj},'kpl')
                parmst(parmi)=parms(parmi)*(1+ef(kk));
            else
                fdvt.knownvals(kpi)=fdvr.knownvals(kpi)*(1+ef(kk));
            end
            [EV, IV, vols]=P2L3s(parmst,fdvt);
            spyr=IV*vols(1) + EV(1,:)*vols(2) + EV(2,:)*vols(3);
            slac=EV(3,:)*vols(3);
            aucp= trapz(spyr.*fdvt.TR);
            aucl= trapz(slac.*fdvt.TR);
            aucs(ii,jj,kk)=aucl/aucp;
        end
    end
end

%% Plot
fig3 = figure('Name', 'Figure 3');
t = tiledlayout(1,3, 'TileSpacing','compact');
% fontsize = 7.5;
fontsize = 15;

fig3.Units = 'inches';
fig3.Position = [2.4583, 3.6771, 10.9583, 5.0208];

red = [0.735, 0.078, 0.184];
green = [0.194, 0.529, 0.025];
yellow = [0.929, 0.694, 0.125]; 

labels = {'A', 'B', 'C'};
y_lower = 100;
y_upper = 0; 
for ii=1:length(testparms)
    ax = nexttile;
    avp = squeeze(aucs(:,ii,:));
    %normalize to value when ef=0:
    for jj=1:3
        avp(jj,:)=avp(jj,:)/avp(jj,(nefpoints+1)/2);
    end

    p1 = plot(ef*100, avp(1,:),'-o', 'color',red,...
            'LineWidth',2, 'DisplayName','Pyr_{iv} Dominant');
    hold on;
    p2 = plot(ef*100, avp(3,:),'-o', 'color',green,...
            'LineWidth',2, 'DisplayName','Pyr_{c}  Dominant',...
            'MarkerFaceColor', green);
    p3 = plot(ef*100, avp(2,:),'-x', 'color',yellow,...
            'LineWidth',2, 'DisplayName','Pyr_{ee} Dominant',...
            'MarkerSize',fontsize*0.6);
    hold off; 
       
    y_lower = min( [y_lower,min(avp,[],'all')] );
    y_upper = max( [y_upper,max(avp,[],'all')] );

    axis tight
    grid on
    
    add_subpanel_label(labels{ii}, ax, fontsize);

    set(gca,'fontsize',fontsize*0.8)
    if strcmp(testparms{ii},'kve')
        xls='k_{ve}';legloc='Northwest';
    elseif strcmp(testparms{ii},'kecp')
        xls='k_{ecP}';legloc='Northwest';
    elseif strcmp(testparms{ii},'kpl')
        xls='k_{PL}';legloc='Northwest';
    elseif strcmp(testparms{ii},'vb')
        xls='v_{b}';legloc='Northeast';
    end        
    xlabel(sprintf('%% Variation in %s',xls),'fontsize',fontsize)
end

% Set all ylims the same for comparison
y_bounds = [y_lower, y_upper];
nexttile(1);
ylim(y_bounds);
nexttile(2);
ylim(y_bounds);
nexttile(3); 
ylim(y_bounds);

nexttile(1);
legend([p1,p3,p2], 'Location',legloc, 'FontSize',fontsize*0.78);
ylabel('Relative Change in AUC Ratio', 'FontSize',fontsize);

% Export as MATLAB fig file for later editing
savefig(fig3, '../Media/figure3.fig');

% Export as vector image for publication
exportgraphics(fig3, '../Media/figure3.eps',...
    'Resolution', 600,...
    'BackgroundColor', 'white',...
    'ContentType', 'vector');

% Export as pixel image for preview
exportgraphics(fig3, '../Media/figure3.png',...
    'Resolution', 600,...
    'BackgroundColor', 'white');


% Local function for adding subpanel labels
function add_subpanel_label(label, axis, fontsize, shift)
    if nargin ~= 3
        xshift = shift(1);
        yshift = shift(2);
    else
        xshift = 0;
        yshift = 0; 
    end

    pos = axis.Position;

    annotation('textbox',[pos(1)+xshift-0.01, pos(2)+pos(4)-0.03+yshift, 0.03, 0.03],...
        'String',label, 'FontSize',fontsize*2, 'FontWeight','bold',...
        'EdgeColor','none', 'VerticalAlignment','bottom');
end

