function [fgfh] = registeredoverlay2(name,fgim,fgh,usfac,bgim,bgh,bgxi,bgyi,anflag,roimask)

% Variables from overlay that we need for registration:
%  Xgh.ImagePositionPatient
%  Xgh.Rows or h.Height?
%  Xgh.Columns or h.Width?
%  Xgh.ReconstructionDiameter

%Set position location for upper left and lower right voxels
% Note: both DICOM and Matlab reference the position of the voxel center.
% for now assume axial
% center of first and last voxel in each row (x) and column (y):

%Are we zooming into a subregion of the registered data?
if ~exist('bgxi','var')
    bgxi=[1 bgh.Columns];
end
if ~exist('bgyi','var')
    bgyi=[1 bgh.Rows];
end
if ~exist('roimask','var') % outline certain voxels
    roimask=zeros(size(fgim));
end

%Are we adding annotations?
if ~exist('anflag','var')
    anflag=0;
end

%Get vector of position for all rows, columns of background image
bgxv=bgh.ImagePositionPatient(1)+(bgh.ReconstructionDiameter*double(0:bgh.Columns-1)/double(bgh.Columns));
bgyv=bgh.ImagePositionPatient(2)+(bgh.ReconstructionDiameter*double(0:bgh.Rows-1)/double(bgh.Rows));
%Position of voxels of zoomed region:
bgx=bgxv(bgxi);
bgy=bgyv(bgyi);

%Same for foreground
fgxv=fgh.ImagePositionPatient(1)+(fgh.ReconstructionDiameter*double(0:fgh.Columns-1)/double(fgh.Columns));
fgyv=fgh.ImagePositionPatient(2)+(fgh.ReconstructionDiameter*double(0:fgh.Rows-1)/double(fgh.Rows));
%Find indices for foreground that fit within background:
fgxi=[1 fgh.Columns];
while fgxv(fgxi(1))<(bgx(1)+fgh.ReconstructionDiameter/double(fgh.Columns)/2/usfac)
    fgxi(1)=fgxi(1)+1;
end
while fgxv(fgxi(2))>(bgx(2)-fgh.ReconstructionDiameter/double(fgh.Columns)/2/usfac)
    fgxi(2)=fgxi(2)-1;
end
fgx=fgxv(fgxi);
fgyi=[1 fgh.Rows];
while fgyv(fgyi(1))<(bgy(1)+fgh.ReconstructionDiameter/double(fgh.Rows)/2/usfac)
    fgyi(1)=fgyi(1)+1;
end
while fgyv(fgyi(2))>(bgy(2)-fgh.ReconstructionDiameter/double(fgh.Rows)/2/usfac)
    fgyi(2)=fgyi(2)-1;
end
fgy=fgyv(fgyi);

%Take care of resampling the overlay here, so we don't get confused on
% registration:
if usfac==1 % Use native resolution of overlay
    fgimint=fgim(fgyi(1):fgyi(2),fgxi(1):fgxi(2));
else
    fgimint=interp2(linspace(fgx(1),fgx(2),fgxi(2)-fgxi(1)+1),... % x positions for native image
        linspace(fgy(1),fgy(2),fgyi(2)-fgyi(1)+1)',... % y positions for native image
        fgim(fgyi(1):fgyi(2),fgxi(1):fgxi(2)),...
        linspace(fgx(1),fgx(2),(fgxi(2)-fgxi(1)+1)*usfac),... % interpolated x positions
        linspace(fgy(1),fgy(2),(fgyi(2)-fgyi(1)+1)*usfac)',... % interpolated y positions
        'spline');
    %Note here that location of corner voxels are the same as the original
    %data
end

%Set alpha for the foreground
fgalpha = fgimint/max(max(fgimint));

%Set up colorbar
cb=zeros(bgyi(2)-bgyi(1)+1,bgxi(2)-bgxi(1)+1);
barthk=floor((bgxi(2)-bgxi(1)+1)/32);
cb(2:end-1,end-barthk:end-1)=repmat(round(linspace(255,1,bgyi(2)-bgyi(1)-1)'),1,barthk);

%Make the magic happen:
bgim2=bgim(bgyi(1):bgyi(2),bgxi(1):bgxi(2));
bgfh=imshow(ind2rgb(floor(255*bgim2/max(max(bgim2))+1),colormap('gray')),'XData',bgx,'YData',bgy);
hold on

fgfh=imshow(ind2rgb(floor(255*fgalpha+1),colormap('hot')),'XData',fgx,'YData',fgy);
set(fgfh,'AlphaData',fgalpha);
cbfh=imshow(ind2rgb(cb,colormap('hot')),'XData',bgx,'YData',bgy);
set(cbfh,'AlphaData',cb>0);
% if(sum(sum(roimask))>0)
%     fgcfh=contour(fgxv,fgyv,roimask,[0.5 0.5]);
% end
hold off

if ~isempty(name)
    text(bgx(1)+0.01*(bgx(2)-bgx(1)),bgy(1)+0.025*(bgy(2)-bgy(1)),...
        'MDACC','color',[1 1 1],'fontsize',14,...
        'fontweight','bold','HorizontalAlignment','Left');
    text(bgx(1)+0.5*(bgx(2)-bgx(1)),bgy(1)+0.025*(bgy(2)-bgy(1)),...
        name,'color',[1 1 1],'fontsize',14,...
        'fontweight','bold','HorizontalAlignment','center');
end

if anflag
mask=fgim>0;
maxfg=max(max(fgim));
text(bgx(1)+0.81*(bgx(2)-bgx(1)),bgy(1)+0.025*(bgy(2)-bgy(1)),...
    sprintf('%5.2e',max(max(fgim))),'color',[1 1 1],'fontsize',12,...
    'fontweight','bold');    
    for jj=fgxi(1):fgxi(2)
        for ii=fgyi(1):fgyi(2)
            if mask(ii,jj)
                if fgim(ii,jj)>0.5*maxfg
                    fontcolor=[0 0 0];
                else
                    fontcolor=[1 1 1];
                end
                if roimask(ii,jj)
                    fontcolor=[0 0 1]; % make key values blue
                    %fontcolor=[0.01 0.4 0.25]; % dark green
                    %fontcolor=[0.0 0.5 0.0]; %
                end
                if maxfg>10
                    text(fgxv(jj),fgyv(ii),sprintf('%4.1f',fgim(ii,jj)),...
                        'Fontsize',14,'FontWeight','bold','Color',fontcolor,...
                        'HorizontalAlignment','center')
                else
                    text(fgxv(jj),fgyv(ii),sprintf('%4.2f',fgim(ii,jj)),...
                        'Fontsize',14,'FontWeight','bold','Color',fontcolor,...
                        'HorizontalAlignment','center')
                end
            end
        end
    end
    
end



