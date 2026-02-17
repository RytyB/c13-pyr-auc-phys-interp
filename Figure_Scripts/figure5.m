close all;
clear; 
clc;

addpath('../Simulation_Functions');

%{
Figure 5

The point: Even with noise, by fitting and calculating the AUC ratio, you
can figure out what someone using a different TE would have measured at a
different site. 

Panel A: noiseless demonstration of the concept.
Panel B: SNR=20 evaluation of robustness to noise with extra fitting parameters
%}

%% User-defined simulation parameters
% n_fits = 1;
n_fits = 500;
n_rand = 50;
snr = 20;
verbose_fitting = 0;
run_sim = 0;


if run_sim
    %% Simulation for Panel A
    
    % Establish a base set of PK parameters
    
    fdv_remote.fitvars = {'kpl', 'VIFScale'};  % representing ground truth data from a remote site
    % Physiological parameters
    fdv_remote.knowns    = {'kve','kecp','vb','vef','R1Lx'};
    fdv_remote.knownvals = [0.0566,0.0293,0.1, 0.5, 0.0393];
    % Physics parameters
    fdv_remote.knowns    = {fdv_remote.knowns{:}, 'T1Pyr', 'T1Lac', 'T2spv','T2spe', 'T2spc', 'T2slc', 'Pe0','Pi0','Li0'};
    fdv_remote.knownvals = [fdv_remote.knownvals,  43,        33,     0.1,   0.055,   0.055,   0.033,    0,    0,    0]; 
    % Parameters for the VIF
    fdv_remote.knowns    = {fdv_remote.knowns{:}, 'Gam1', 'Gam2', 'tdel'}; 
    fdv_remote.knownvals = [fdv_remote.knownvals,  2.8,    4.5,    10]; 
    % Describe acquisition scheme
    fdv_remote.ntp = 64;
    fdv_remote.NSeg = 1;
    fdv_remote.NFlips = fdv_remote.ntp*fdv_remote.NSeg;  % Total number of excitations
    % Describe temporal sampling scheme
    fdv_remote.TE = 0.0219; % in vivo data used 21.9 ms
    fdv_remote.TR = 3 * ones(1,fdv_remote.NFlips);
    fdv_remote.taxis = cumsum(fdv_remote.TR) - fdv_remote.TR(1);
    % Describe the excitation scheme
    fdv_remote.FlipAngle = [20,0; 0,30]*ones(2,fdv_remote.NFlips);
    % Describe the Vascular Input Function (VIF)
    fdv_remote.UseVIF = 1;
    fdv_remote.VIFP = gampdf(fdv_remote.taxis-6, 2.8, 4.5);
    fdv_remote.Name = 'Testing effects of TE/T2s P2L3s';
    % Describe slice profile
    fdv_remote.SPBins = 1;
    fdv_remote.SPWeights = 1;
    % Miscellaneous
    fdv_remote.verbose = 0;
    
    % Fitting parameters
    %    kpl,  VIF ampl
    lb = [0,    0];
    ub = [Inf, Inf];
    jbopts = optimset('display','off');
    
    parms = [0.5371, 1000];
    
    fprintf('Noiseless simulation...\n');
    TEs = (0:5:40) / 1000; % ms
    reported_aucs = zeros(1,length(TEs));
    calcd_aucs = zeros(1,length(TEs));
    corrected_aucs = zeros(1,length(TEs));
    for TEi = 1:length(TEs)
        fdv_remote.TE = TEs(TEi);
        fprintf('\tStarting point %d of %d...\n', TEi, length(TEs));
    
        % Simulate acquiring same physiological data with each TE
        [EV, IV, vols] = P2L3s(parms, fdv_remote);
        fdv_remote.data = [IV*vols(1) + EV(1,:)*vols(2) + EV(2,:)*vols(3);
                    EV(3,:)*vols(3)];
    
        % Fitting loop
        bestresid = Inf;
        bestfits = zeros(1,length(parms));
        for kk = 1:n_fits
            guess = parms*10.*rand(1,2);
            [fits, resid] = lsqnonlin(@(x) P2L3sErr(x,fdv_remote), guess,lb,ub,jbopts);
            if resid < bestresid
                bestresid = resid;
                bestfits = fits;
            end
        end 
    
        [AUCr_calcd,~,~,~,AUCr_remote] = P2L3sAUC(bestfits, fdv_remote, fdv_remote);
        fdv_remote.TE = 0.020;  
        AUCr_corrected = P2L3sAUC(bestfits, fdv_remote, fdv_remote);
    
        reported_aucs(TEi) = AUCr_remote;
        calcd_aucs(TEi) = AUCr_calcd;
        corrected_aucs(TEi) = AUCr_corrected;
    
    end
    
    save('fig5a_data_full_sim.mat',...
        'reported_aucs', 'TEs', 'corrected_aucs', 'reported_aucs', 'calcd_aucs');
    
    
    %% Simulation for Panel B
    fprintf('\nFull simulation...\n');
    TEs = (0:5:40) / 1000; % ms
    reported_aucs = zeros(n_rand,length(TEs));
    calcd_aucs = zeros(n_rand,length(TEs));
    corrected_aucs = zeros(n_rand,length(TEs));
    for TEi = 1:length(TEs)
        fdv_remote.TE = TEs(TEi);
        fprintf('\tStarting point %d of %d...\n', TEi, length(TEs));
    
        [EV, IV, vols] = P2L3s(parms, fdv_remote);
        ground_truth = [IV*vols(1) + EV(1,:)*vols(2) + EV(2,:)*vols(3);
                        EV(3,:)*vols(3)];
    
        for r = 1:n_rand
            % Simulate acquiring same physiological data with each TE
            mu = 0;
            sigma = max(ground_truth, [], 'all') / snr; 
            noise_array = normrnd(mu, sigma, size(ground_truth)); 
            fdv_remote.data = ground_truth + noise_array;
    
            [bestfits,fit_curves,fdv_remote] = fit_3PCs(fdv_remote, n_fits, verbose_fitting);
    
            [AUCr_calcd,~,~,~,AUCr_remote] = P2L3sAUC(bestfits, fdv_remote, fdv_remote);
            fdv_corrected = fdv_remote;
            fdv_corrected.TE = 0.0219;  
            AUCr_corrected = P2L3sAUC(bestfits, fdv_corrected, fdv_corrected);
    
            reported_aucs(r,TEi) = AUCr_remote;
            calcd_aucs(r,TEi) = AUCr_calcd;
            corrected_aucs(r,TEi) = AUCr_corrected;
        end
    
        % figure('Name','Debug');
        % plot(ground_truth(1,:), 'k-');
        % hold on;
        % plot(fdv_remote.data(1,:), 'kx');
        % plot(fit_curves(1,:), 'g-');
        % hold off;
        % drawnow;
    
    end
    
    save('fig5b_data_full_sim.mat',...
        'TEs', 'corrected_aucs', 'reported_aucs', 'calcd_aucs');
end

%% Plotting

fig5 = figure('Name', 'figure5');
tiledlayout(1,1, 'TileSpacing','compact');
fontsize = 14;
markersize = fontsize * 0.8;

fig5.Units = 'inches';
fig5.Position = [7.3250, 2.0083, 5.8333, 3.4500];

% data5a = load('fig5a_data_full_sim.mat');
% TEs = data5a.TEs * 1000; % in ms for display
% ax1 = nexttile; 
% plot(TEs, data5a.reported_aucs, 'kx',...
%     'MarkerSize',markersize, 'LineWidth',1.25, 'DisplayName', 'Integrated AUC Ratio');
% hold on;
% plot(TEs, data5a.calcd_aucs, 'r-',...
%     'LineWidth',1.25, 'DisplayName','Model-Derived AUC Ratio');
% plot(TEs, data5a.corrected_aucs, 'b:',...
%     'LineWidth',2, 'DisplayName','TE Corrected AUC Ratio');
% hold off;
% set(gca, 'FontSize',fontsize*0.8);
% ylabel('AUC Ratio', 'FontSize',fontsize, 'FontWeight','bold');
% % legend('FontSize',fontsize*0.8, 'Location','southwest');
% grid on;
% % ylim([0.15, max(data5a.reported_aucs)]);

data5b = load('fig5b_data_full_sim.mat');
TEs = data5b.TEs;
calcd_aucs = data5b.calcd_aucs;
reported_aucs = data5b.reported_aucs;
corrected_aucs = data5b.corrected_aucs;
ax2 = nexttile; 
plot(TEs*1000, mean(data5b.reported_aucs, 1), 'kx',...
    'MarkerSize',markersize, 'LineWidth',1.25, 'DisplayName', 'Integrated AUC Ratio');
hold on;
errorbar(TEs*1000, mean(data5b.calcd_aucs,1), std(data5b.calcd_aucs,[],1), 'r-',...
    'LineWidth',1.25, 'DisplayName','Model-Derived AUC Ratio');
errorbar(TEs*1000, mean(data5b.corrected_aucs,1),std(data5b.corrected_aucs,[],1), 'b:',...
    'LineWidth',2, 'DisplayName','TE Corrected AUC Ratio');
hold off;
set(gca, 'FontSize',fontsize*0.8);
ylabel('AUC Ratio', 'FontSize',fontsize, 'FontWeight','bold');
xlabel('Echo Time (ms)', 'FontSize',fontsize, 'FontWeight','bold');
legend('FontSize',fontsize*0.8, 'Location','southwest');
grid on;



%% Local function for fitting
function [popt, fit_curves, fdv_modded, nmse, convg_study] = fit_3PCs(fdv, n_fits, verbose)
    jbopts = optimset('display', 'off'); 
    best_resid = Inf;
    best_fits = 0;
    
    if verbose
        fprintf('Starting fitting with %d repetitions...', n_fits);
        tic;
    end

    % % For fitting T1s without extra loss factor
    % fdv.fitvars   = {'kpl', 'VIFScale', 'kve', 'kecp', 'vb', 'vef', 'T1Pyr', 'T1Lac'};
    % fdv.knowns    = {'klp', 'kecl', 'R1Lx', 'T2spv', 'T2spe', 'T2spc', 'T2slc', 'Pe0', 'Pi0', 'Le0', 'Li0', 'FAp', 'FAl'};
    % fdv.knownvals = [ 0.0,   0.0,    0.0,     0.1,    0.055,   0.055,   0.033,    0,     0,     0,     0,     20,   30];
    % UB            = [ 10,      Inf,       10,    10,    1,     1,     Inf,      Inf ];
    % LB            = [  0,       0,         0,     0,    0,     0,      0,        0 ];
    % p0            = [ 0.5,     1000,     0.07,  0.02,   0.1,  0.5,    43,       33 ];
    
    % % For fitting extra loss factor without T1s
    fdv.fitvars   = {'kpl', 'VIFScale', 'kve', 'kecp', 'R1Lx'};
    UB            = [  1,      Inf,       1,     1,       1 ];
    LB            = [  0,       0,         0,    0,       0 ];
    p0            = [ 0.5,     1000,     0.07,  0.02,     0 ];

    convergence = 0;
    resids = Inf; 
    for ii=1:n_fits
        guess = (p0*9).*rand(1,length(p0));
        [fits, resid] = lsqnonlin(@(x) P2L3sErr(x, fdv), guess, LB, UB, jbopts);

        if resid < best_resid
            best_resid = resid;
            best_fits = fits;
            convergence = [convergence, ii];
            resids = [resids, best_resid];
        end
    end
    convg_study = [convergence; resids];

    popt = best_fits;
    [EV, IV, vols] = P2L3s(popt, fdv);
    fit_curves = [
        IV*vols(1) + EV(1,:)*vols(2) + EV(2,:)*vols(3);  % Pyruvate signal
        EV(3,:)*vols(3)                                  % Lactate signal
    ]; 
    fdv_modded = fdv;
    nmse = calc_nmse(fdv.data, fit_curves);

    if verbose
        fprintf('finished. Elapsed time: %.2f min\n', toc()/60);
    end
end