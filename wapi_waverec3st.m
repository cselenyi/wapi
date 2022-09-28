function a = wapi_waverec3st(wd3Fname,s,varargin)
%wapi_waverec3st        Multi-level 3D wavelet reconstruction
%
% a = wapi_waverec3st(c,s,rh,rg)
% a = wapi_waverec3st(c,s,rh,rg,rhz,rgz)
% a = wapi_waverec3st(wd3Fname,[],rh,rg)
% a = wapi_waverec3st(wd3Fname,[],rh,rg,rhz,rgz)
%
% wapi_waverec3st performs a multi-level 3D wavelet reconstruction using
% specific reconstruction filters. Both dyadic and stationary
% reconstruction is supported. Presence of stationary wavelet decomposition
% structure is detected automatically.
%
% Inputs:
%
% c             The decomposition coefficients, see <a href="matlab: help wapi_wavedec3st">wapi_wavedec3st</a> for 
%               details.
%
% s             The decomposition bookkeeping matrix, see <a href="matlab: help wapi_wavedec3st">wapi_wavedec3st</a>
%               for details.
%
% wd3Fname      If the decomposition data is in a .wd3 file then you can
%               provide the name of it as the first input. In this case the
%               second argument must be an empty matrix
%
% rh            Reconstruction low-pass filter.
%
% rg            Reconstruction high-pass filter.
%
% rhz           Optional. If anisotrop WT decomposition was used then a
%               separate reconstruction low-pass filter must be specified
%               for the Z-axis (X & Y will always use filter rh and rg). If
%               rhz is specified then rgz must be also specified.
%
% rgz           Optional. Reconstruction high-pass filter to be used along
%               Z-axis (see rhz above).
%
%
% Output:
%
% a             Reconstructed 3D array (image).
%

%
%     Copyright (C) 2022 by Zsolt Cselényi
%
%     The WAPI toolbox (collection of functions listed under the heading
%     "Proper WAPI toolbox functions (toolbox manifest)" in wapi_help.m
%     file) is free software: you can redistribute it and/or modify it
%     under the terms of the GNU General Public License as published by the
%     Free Software Foundation, either version 3 of the License, or (at
%     your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%     Author: Zsolt Cselényi
%     e-mail: zsolt.cselenyi@ki.se
%
%     Version 2022-04-14

if ~any((3:7)==nargin)
    error('Wrong number of inputs');
end
if ~any((0:1)==nargout)
    error('Wrong number of outputs');
end

if isempty(s)
    [offs,numfr,lc,s]=wapi_readwd3(wd3Fname); %#ok<*ASGLU>
else
    lc=sum(prod([s(1:end-1,:) [1;repmat(7,size(s,1)-2,1)]],2)); %#ok<NASGU>
end

rmax = size(s,1);
nmax = rmax-2;
anisotrop_filter=0;
if ischar(varargin{1})
    [LoF_R1,HiF_R1] = wfilters(varargin{1},'r'); % must have Wavelet Toolbox for this
    if(nargin>3 && ischar(varargin{2}))
        [LoF_R2,HiF_R2] = wfilters(varargin{2},'r'); 
        anisotrop_filter=1;
    else
        LoF_R2=LoF_R1;
        HiF_R2=HiF_R1;
    end
else
    LoF_R1 = varargin{1}; HiF_R1 = varargin{2}; 
    if(nargin == 4) || (nargin>4 && ischar(varargin{3}))
        LoF_R2=LoF_R1;
        HiF_R2=HiF_R1;
    else
        LoF_R2 = varargin{3}; HiF_R2 = varargin{4}; 
        anisotrop_filter=1;
    end
end

st=0;
if ischar(varargin{end}) 
    if strcmp(varargin{end},'st')
        st=1;
    end
else
    if all(sum(diff(s))==0)
        disp('Stationary decomposition structure detected.');
        st=1; % trimmed stationary transform with trimmed s matrix
    end
end

n=0; % limit reconstruction to reach this level+1 only from deepest

[fst,lst]=wapi_idxcoef3('app',s,1);
sz=s(1,:);
if ischar(wd3Fname)
    a=reshape(wapi_readwd3(wd3Fname,fst,lst,1),sz);
else
    a=reshape(wd3Fname(fst:lst),sz);
end

if st
    to=nmax;

    for p=2:to
        LoF_R2=dyadup(LoF_R2);
        LoF_R1=dyadup(LoF_R1);
        HiF_R2=dyadup(HiF_R2);
        HiF_R1=dyadup(HiF_R1);
    end
end
  
rm   = rmax+1;
fillArrayFromFile=@wapi_mexFillArrayFromFile;
if version('-release')<14
    error('Pre-release 14 is not supported.');
else
    idwt3mx=@wapi_idwt3mex_v7;
    idwt3zmx=@wapi_idwt3zmex_v7;
end
idwt3z_xlow_mx=@wapi_idwt3z_xlow_mex;
idwt3z_xhigh_mx=@wapi_idwt3z_xhigh_mex;
split_x_low_high=true;

if ~st
    split_x_low_high=false;
end

[nm,mx,endian]=computer;
if strcmp(endian,'B')
    byteSwap=1;
else
    byteSwap=0;
end
cnames_xlow={'app','h_l','a_h','h_h'};
cnames_xhigh={'v_l','d_l','v_h','d_h'};
for p=nmax:-1:n+1
    if split_x_low_high
        sz=[s(end-p,:) 4];
        c=zeros(sz);
        c(:,:,:,1)=a;
        a(:)=0; %=zeros(sz(1:3));
        for i=2:4
            [fst,lst]=wapi_idxcoef3(cnames_xlow{i},s,p);
            offsetArr=prod(sz(1:3))*(i-1)+1;
            fillLength=lst-fst+1;
            offsetFile=fst*8; % no -1 since first double must be skipped
            if ischar(wd3Fname)
                fillArrayFromFile(c,offsetArr,fillLength,[wd3Fname '1'],offsetFile,byteSwap);
            else
                c(offsetArr:offsetArr+fillLength-1)=wd3Fname(fst:lst);
            end
        end
        disp('Calling idwt3z_xlow_mx');
        drawnow;
        if st
            if p<nmax
                LoF_R2=dyaddown(LoF_R2);
                LoF_R1=dyaddown(LoF_R1);
                HiF_R2=dyaddown(HiF_R2);
                HiF_R1=dyaddown(HiF_R1);
            end
            if anisotrop_filter
                a = feval(idwt3z_xlow_mx,c,LoF_R1,HiF_R1,LoF_R2,HiF_R2,a,'st');
            else
                a = feval(idwt3z_xlow_mx,c,LoF_R1,HiF_R1,LoF_R1,HiF_R1,a,'st');
            end
        else
            if anisotrop_filter
                a = feval(idwt3z_xlow_mx,c,LoF_R1,HiF_R1,LoF_R2,HiF_R2,a);
            else
                a = feval(idwt3z_xlow_mx,c,LoF_R1,HiF_R1,LoF_R1,HiF_R1,a);
            end
        end
        for i=1:4
            [fst,lst]=wapi_idxcoef3(cnames_xhigh{i},s,p);
            offsetArr=prod(sz(1:3))*(i-1)+1;
            fillLength=lst-fst+1;
            offsetFile=fst*8; % no -1 since first double must be skipped
            if ischar(wd3Fname)
                fillArrayFromFile(c,offsetArr,fillLength,[wd3Fname '1'],offsetFile,byteSwap);
            else
                c(offsetArr:offsetArr+fillLength-1)=wd3Fname(fst:lst);
            end
        end
        disp('Calling idwt3z_xhigh_mx');
        drawnow;
        if st
            if anisotrop_filter
                a = feval(idwt3z_xhigh_mx,c,LoF_R1,HiF_R1,LoF_R2,HiF_R2,a,'st');
            else
                a = feval(idwt3z_xhigh_mx,c,LoF_R1,HiF_R1,LoF_R1,HiF_R1,a,'st');
            end
        else
            if anisotrop_filter
                a = feval(idwt3z_xhigh_mx,c,LoF_R1,HiF_R1,LoF_R2,HiF_R2,a);
            else
                a = feval(idwt3z_xhigh_mx,c,LoF_R1,HiF_R1,LoF_R1,HiF_R1,a);
            end
        end
        clear c
    else
        [fst,lst]=wapi_idxcoef3('c',s,p);
        if 1
            sz=[s(end-p,:) 8];
            tmp=a;
            a=zeros(sz);
            a(:,:,:,1)=tmp;
            clear tmp
            offsetArr=prod(sz(1:3))+1;
            fillLength=lst-fst+1;
            offsetFile=fst*8;
            if ischar(wd3Fname)
                fillArrayFromFile(a,offsetArr,fillLength,[wd3Fname '1'],offsetFile,byteSwap);
            else
                a(offsetArr:offsetArr+fillLength-1)=wd3Fname(fst:lst);
            end
        else
            sz=[s(end-p,:) 7]; %#ok<UNRCH>
            if ischar(wd3Fname)
                tmp=reshape(wapi_readwd3(wd3Fname,fst,lst,1),sz);
            else
                tmp=reshape(wd3Fname(fst:lst),sz);
            end
            a(:,:,:,2:8)=tmp;
        end
        if st
            if p<nmax
                LoF_R2=dyaddown(LoF_R2);
                LoF_R1=dyaddown(LoF_R1);
                HiF_R2=dyaddown(HiF_R2);
                HiF_R1=dyaddown(HiF_R1);
            end
            if anisotrop_filter
                a = feval(idwt3zmx,a,LoF_R1,HiF_R1,LoF_R2,HiF_R2,s(rm-p,:),'st');
            else
                a = feval(idwt3mx,a,LoF_R1,HiF_R1,s(rm-p,:),'st');
            end
        else
            if anisotrop_filter
                a = feval(idwt3zmx,a,LoF_R1,HiF_R1,LoF_R2,HiF_R2,s(rm-p,:));
            else
                a = feval(idwt3mx,a,LoF_R1,HiF_R1,s(rm-p,:));
            end
        end
    end
    fprintf(1,'Level %d done.\n',p);
end

function y=dyadup(x)

l = 2*length(x)+1;
y = zeros(1,l);
y(2:2:l) = x;

function y=dyaddown(x)

y = x(2:2:end);

