%% Determine baseline AUC Ratio and Pyruvate Volume Fraction in HNSCC patients:
clear 

addpath('../Simulation_Functions');
dsets={'../Data/Subj01a','../Data/Subj01b','../Data/Subj02','../Data/Subj03','../Data/Subj04'};

%{
Figure 7

The Point: Here is the novel analysis applied to several human datasets.

Note: This script will spit out a bunch of PNGs that need to be stitched
together to make the version of Figure 7 that appears in the manuscript. 
%}


%% Analyze data, generate PVF mask:

voxelcount=0;
PyrCompartment=[0 0 0];

%p2u=[1 1 2 3 7];
%For collage:
%r2u=[1 2 1 1 1]; %rep
s2u=[5 5 4 4 5]; %slice
%figure(8);clf
%tlh=tiledlayout(2,length(p2u),'TileSpacing','Compact','Padding','Compact');

zi=37;
auc_m=[];
auc_c=[];
pvf_c=[];

for ii=1:length(dsets)
    load(dsets{ii})
    [NPE NOP NSL]=size(tmpmask);
    PyrCompMask=zeros(NPE,NOP,NSL,3);
    aucc=zeros(NPE,NOP,NSL);
    aucm=zeros(NPE,NOP,NSL);
    for kk=1:NPE
        for mm=1:NOP
            for nn=1:NSL
                if tmpmask(kk,mm,nn)
                    voxelcount=voxelcount+1;
                    fdv3sV2.knownvals(kpli)=kplcs2s(kk,mm,nn,5);
                    fitparms=squeeze(kplcs2parms(kk,mm,nn,5,:));
                    % [~,aucc_tmp,pvf,~]=P2L3sAUC(fitparms,fdv3sV2);
                    [aucc_tmp,pvf] = P2L3sAUC(fitparms, fdv3sV2);
                    auc_c=[auc_c aucc_tmp];
                    aucc(kk,mm,nn)=aucc_tmp;
                    [~,pvfi]=max(pvf);
                    pvf_c=[pvf_c pvfi];
                    aucm(kk,mm,nn)=trapz(slac(kk,mm,nn,laci(1):laci(2)))/trapz(spyr(kk,mm,nn,pyri(1):pyri(2)));
                    auc_m=[auc_m aucm(kk,mm,nn)];
                    PyrCompartment(pvfi)=PyrCompartment(pvfi)+1;
                    PyrCompMask(kk,mm,nn,pvfi)=1;
                    % if (ii==1)&&(jj==1)&&(kk==8)&&(mm==11)&&(nn==5)
                    %     fprintf('First Rep PVF=(%5.3f,%5.3f,%5.3f), Pat %s TP %d Rep %d Vox (%d,%d,%d)!\n',...
                    %         pvf,pat(pn).Name,tp,jj,kk,mm,nn);
                    % end
                    % if pvfi==1
                    %     fprintf('IV Dominant PVF=(%5.3f,%5.3f,%5.3f), Pat %s TP %d Rep %d Vox (%d,%d,%d)!\n',...
                    %         pvf,pat(pn).Name,tp,jj,kk,mm,nn);
                    % elseif pvfi==3
                    %     fprintf('IC Dominant PVF=(%5.3f,%5.3f,%5.3f), Pat %s TP %d Rep %d Vox (%d,%d,%d)!\n',...
                    %         pvf,pat(pn).Name,tp,jj,kk,mm,nn);
                    %end
                end %tmpmask
            end %NSL
        end %NOP
    end %NPE
    %if jj==r2u(ii) % generate the overlay:
    fgim=squeeze(aucc(:,:,s2u(ii)));
    f71h=figure(71);%(ii) %(1,ii)
    registeredoverlay2('',fgim,hph,1,bgim,bgh,[37 476],[37 476],1);
    f71h.Position=[1950 -250 940 820];
    exportgraphics(gca,sprintf('Fig7_AUCc_%d.png',ii));
    %
    fgim=squeeze(aucm(:,:,s2u(ii)));
    f72h=figure(72);%(ii) %(1,ii)
    registeredoverlay2('',fgim,hph,1,bgim,bgh,[37 476],[37 476],1);
    f72h.Position=[2900 -250 940 820];
    exportgraphics(gca,sprintf('Fig7_AUCi_%d.png',ii));
    %
    f73h=figure(73); %nexttile(ii+5) %(2,ii)
    imshow(bgim(zi:end-zi+1,zi:end-zi+1),[],'XData',bgx(zi:end-zi+1),'Ydata',bgy(zi:end-zi+1));
    hold on
    contour(bgx,bgy,tumroi,[0.5 0.5],'c-','linewidth',1);
    % % %contour(hpx,hpy,PyrCompMask(:,:,s2u(ii),2),[0.5 0.5],'y-','linewidth',1);
    % % xedge=[hpx-(hpx(2)-hpx(1))/2 hpx(end)+(hpx(2)-hpx(1))/2];
    % % yedge=[hpy-(hpy(2)-hpy(1))/2 hpy(end)+(hpy(2)-hpy(1))/2];
    % % [xc,yc]=meshgrid(xedge,yedge);
    % % contour(xc,yc,PyrCompMask(:,:,s2u(ii),2),[0.5 0.5],'y-','linewidth',1);
    pcme=squeeze(PyrCompMask(:,:,s2u(ii),2));
    % pcme2=interp2(hpx,hpy',pcme,bgx,bgy');
    % pcme2(isnan(pcme2))=0;
    % contour(bgx,bgy,pcme2,[0.5 0.5],'y-','linewidth',1);
    hrpcm=zeros([size(bgim) 3]);
    for yi=1:size(bgim,1);
        [~,hpyi]=min(abs(bgy(yi)-hpy));
        for xi=1:size(bgim,2)
            [~,hpxi]=min(abs(bgx(xi)-hpx));
            hrpcm(yi,xi,:)=PyrCompMask(hpyi,hpxi,s2u(ii),:);
        end
    end
    contour(bgx,bgy,squeeze(hrpcm(:,:,2)),[0.5 0.5],'y-','linewidth',2);
    if sum(sum(squeeze(hrpcm(:,:,1))))>0
        contour(bgx,bgy,squeeze(hrpcm(:,:,1)),[0.5 0.5],'r-','linewidth',2);
    end                    
    if sum(sum(squeeze(hrpcm(:,:,3))))>0
        contour(bgx,bgy,squeeze(hrpcm(:,:,3)),[0.5 0.5],'g-','linewidth',2);
    end
    hold off
    f73h.Position=[3850 -250 940 820];
    exportgraphics(gca,sprintf('Fig7_PVCROI_%d.png',ii))
    fprintf('Press any key to continue...\n'); 
    pause;
    fprintf('\n');
    %end
end %

PyrCompartment/voxelcount
% [2,54,2]/58=[0.0345 0.9310 0.0345]

% Need to include 006 slice 5 for IV dominant, 
%                 012 slice 4-5 for IC dominant.

%exportgraphics(gcf,'BaselinePVF_HNSCC.png')

%% regression

cc=corrcoef(auc_m,auc_c);
Correlation_Coefficient=cc(1,2)

lm=fitlm(auc_m,auc_c)
%lm.Coefficients.pValue
R2Adj=lm.Rsquared.Adjusted
m=lm.Coefficients.Estimate 
if length(m)>1
    b=m(1);
    m=m(2);
else
    b=0;
end

figure(8)
plot(auc_m,auc_c,'kx',[0 2.1],b+[0 2.1]*m,'k--','linewidth',2,'markersize',10);
hold on
plot(auc_m(pvf_c==1),auc_c(pvf_c==1),'ro',...
    auc_m(pvf_c==3),auc_c(pvf_c==3),'go','linewidth',1,'markersize',10);
hold off
ylabel('AUC Calculated from Model Parms')
xlabel('AUC Integrated over Signal Curves')
grid on
axis([0 2.1 0 2.1]);
refline(1,0)
exportgraphics(gca,'Fig8aV1.png');
save('Fig8Data.mat','auc_c','auc_m','lm','Correlation_Coefficient')
