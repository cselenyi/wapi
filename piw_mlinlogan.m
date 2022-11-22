function [B,Bstd] = piw_mlinlogan(V,time,ifunc,nframes,w,linreg,k2r, aucthr)
%piw_mlinlogan          Multi-linear (reference-)Logan plot
%
% [B,Bstd] = piw_mlinlogan(V,time,ifunc,nframes,w,[linreg],[k2r],[aucthr])
%
% Calculates linear or multi-linear (reference-)Logan fit. It uses fast
% multithreaded code (piw_mlinlogan_mex) for performing the calculations.
%
% Inputs:
%
% V             kinetic curves (num_frames x num_voxels matrix)
%
% ifunc         plasma or reference region input fuction (2 columns for
%               time and radioactivity concentration)
%
% time          timing of the frames
%
% nframes		number of frames to be used for the linear regression
%
% w             weights to be used for the linear regression
%
% linreg		flag for Linear Regression
%               linreg== 0 (default): Multilinear Regression
%               linreg~= 0: Standard Linear Regression (Logan plot)
%
% k2r           rate constant for Reference Logan.
%               >0		correction applied
%               Otherwise	no correction (default)
%
% aucthr        Optional. If given it specifies the AUC threshold to be
%               used for getting valid coefficient TAC's for Logan
%               calculation (result for below threshold coeff's are set to
%               0 (std=1). Default is no thresholding.
%
% Outputs:
%
% B             Logan slope = volume of distribution (ratio) 
%
% Bstd          Corresponding standard error of the slope coefficients
%

%
% piw_mlinlogan:
%
% Based on PIW_LOGAN from PIWAVE package.
% Modifed by Zsolt Cselényi. 2022-11-04
%
% piw_mlinlogan_mex:
% Derived from PIW_LOGAN from PIWAVE package.
% Written by Zsolt Cselényi. 2022-04-14
%
% PIWave is copyright, distributed under the GNU General Public Licence.
% Please see piw_licence.man for details.
% 
% PIWave written by 
% Federico Turkheimer
%

[d,numVox]=size(V);
if nargin < 6, linreg=0; end
if nargin < 7, k2r=0; end
if nargin < 8, aucthr=[]; end

%	Defining the input function
act_if	= ifunc(:,2);
time_if	= ifunc(:,1);
act_if	= [0;act_if];
time_if= [0;time_if];
time= [0;time];
act_if_int=zeros(size(act_if));
for k=2:length(act_if)
    act_if_int(k)=trapz(time_if(1:k),act_if(1:k));
end;
% act_if_int(1)=0;

if (k2r>0)
   act_if_int = act_if_int+(act_if/k2r);
end
		
% 	Interpolation of the input function
if isequal(time_if,time)
    act_if_i=act_if;
    act_if_int_i=act_if_int;
else
    if exist('octave_config_info','builtin')==5
        act_if_i	= interp1(time_if,act_if,time);
        act_if_int_i= interp1(time_if,act_if_int,time);
    else
        act_if_i	= interp1q(time_if,act_if,time);
        act_if_int_i= interp1q(time_if,act_if_int,time);
    end
end

B=zeros(1,numVox);
Bstd=ones(1,numVox);

if ~isempty(aucthr)
    if numel(aucthr)>1
        idx=find(aucthr);
    else
        auc=sum(V,1); % trapz(time(2:end),V);
        idx = find(abs(auc)>=aucthr);
        fprintf(1,'Percent of valid: %.2f\n',length(idx)./numVox.*100);
    end
    if isempty(idx)
        return;
    end
else
    idx = 1:numVox;
end

V1(2:d+1,:)=V(:,idx);

[B(1,idx),Bstd(1,idx)]=piw_mlinlogan_mex(V1,time,act_if_i,act_if_int_i,nframes,w,linreg);

doPause=0;
bidx=isnan(B(1,idx));
if any(bidx)
    bidx2=sum(V1(end-nframes+1:end,bidx).*repmat(w(end-nframes+1:end,1),[1 sum(double(bidx))]),1)==0;
    tmp=false(size(bidx));
    tmp(bidx)=bidx2;
    bidx2=tmp;
    clear tmp
    B(1,idx(bidx2))=0;
    Bstd(1,idx(bidx2))=1;
    bidx=find(bidx & ~bidx2);
    for i=1:length(bidx)
        B(1,idx(bidx(i)))=0;
        Bstd(1,idx(bidx(i)))=1;
        disp('NaNs in B!');
        drawnow
        disp('V1')
        V1(:,bidx(i))
        disp('act_if_i');
        act_if_i %#ok<NOPRT>
        disp('act_if_int_i');
        act_if_int_i %#ok<NOPRT>
        if doPause
            pause
            doPause=0;
        end
    end
end