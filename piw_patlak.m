function [B,Bstd,A] = piw_patlak(V,time,plasma,nframes,w,aucthr)
%piw_patlak             Patlak plot
%
% [B,Bstd,A] = piw_patlak(V,time,ifunc,nframes,w,[aucthr])
%
% Calculates Patlak fit. It uses fast multithreaded code (piw_patlak_mex) for performing the calculations.
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
% aucthr        Optional. If given it specifies the AUC threshold to be
%               used for getting valid coefficient TAC's for Logan
%               calculation (result for below threshold coeff's are set to
%               0 (std=1). Default is no thresholding.
%
% Outputs:
%
% B             Patlak slope = Ki (ratio)
%
% Bstd          Corresponding standard error of the slope coefficients
%
% A             Patlak intercept = V_ND volume of non-displaceable distribution (ratio)

%
% piw_patlak:
%
% Based on PIW_PATLAK from PIWAVE package.
% Modifed by Zsolt Cselényi. 2022-11-04
%
% piw_patlak_mex:
% Derived from PIW_PATLAK from PIWAVE package.
% Written by Zsolt Cselényi. 2022-11-04
%
% PIWave is copyright, distributed under the GNU General Public Licence.
% Please see piw_licence.man for details.
% 
% PIWave written by 
% Federico Turkheimer
%

[d,numVox]=size(V);
if nargin < 6, aucthr=[]; end

%	Defining the input function
glu	= plasma(:,2);
timep	= plasma(:,1);
glu	= [0;glu];
timep= [0;timep]; % should be in minutes
time= [0;time]; % should be in minutes
gluint=zeros(size(glu));
for k=2:length(glu)
    gluint(k)=trapz(timep(1:k),glu(1:k));
end;
% gluint(1)=0;
		
% 	Interpolation of the input function
if isequal(timep,time)
    glu_i=glu;
    gluint_i=gluint;
else
    if exist('octave_config_info','builtin')==5
        glu_i	= interp1(timep,glu,time);
        gluint_i= interp1(timep,gluint,time);
    else
        glu_i	= interp1q(timep,glu,time);
        gluint_i= interp1q(timep,gluint,time);
    end
end

B=zeros(1,numVox);
Bstd=ones(1,numVox);
if nargout==3
	A=zeros(1,numVox);
end

if ~isempty(aucthr)
    if numel(aucthr)>1
        idx=find(aucthr);
    else
        if exist('octave_config_info','builtin')==5
            auc=trapz(repmat(time(2:end),[1 numVox]),V);
        else
            auc=trapz(time(2:end),V); % dunno why it was replaced with this before: sum(V,1); ???
        end
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

if exist('octave_config_info','builtin')==5
    [B(1,idx),Bstd(1,idx)]=piw_patlak_oct(V1,time,glu_i,gluint_i,nframes,w);
else
	if nargout==3
	    [B(1,idx),Bstd(1,idx),A(:,idx)]=piw_patlak_mex(V1,time,glu_i,gluint_i,nframes,w);
	else
	    [B(1,idx),Bstd(1,idx)]=piw_patlak_mex(V1,time,glu_i,gluint_i,nframes,w);
	end
end

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
	if nargout==3
		A(1,idx(bidx2))=0;
	end
    bidx=find(bidx & ~bidx2);
    for i=1:length(bidx)
        B(1,idx(bidx(i)))=0;
        Bstd(1,idx(bidx(i)))=1;
		if nargout==3
			A(1,idx(bidx2))=0;
		end
        disp('NaNs in B!');
        drawnow
        disp('V1')
        V1(:,bidx(i))
        disp('glu_i');
        glu_i %#ok<NOPRT>
        disp('gluint_i');
        gluint_i %#ok<NOPRT>
        if doPause
            pause
            doPause=0;
        end
    end
end