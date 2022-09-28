function [c,s] = wapi_wavedec3st(x,n,varargin)
%wapi_wavedec3st        Multi-level 3D wavelet decomposition
%
% [c,s] = wapi_wavedec3st(x,n,h,g)
% [c,s] = wapi_wavedec3st(x,n,h,g,st)
% [c,s] = wapi_wavedec3st(x,n,h,g,hz,gz)
% [c,s] = wapi_wavedec3st(x,n,h,g,hz,gz,st)
% [c,s] = wapi_wavedec3st(x,n,...,'save',outfname,...)
%
% [c,s] = wapi_wavedec3st(x,n,..) returns the wavelet decomposition
% of the 3D array x at maximum depth n. Outputs are the decomposition
% vector c and the corresponding bookkeeping matrix s. Both dyadic and
% stationary transform is supported. The filtering step in the wavelet
% decomposition process is performed in the following way for the dyadic
% transform:
%   1. Convolution using the wavelet decomposition filter (low- or
%   high-pass depending axis and subband, see below). The convolution will
%   produce full convolution length, i.e. a data length X and kernel length
%   Y will produce a filtered data of length X+Y-1.
%   2. The filtered data is decimated in the dyadic step: every second
%   sample is thrown away.
% In contrast, the stationary transform works in this way:
%   1. Convolution using the wavelet decomposition filter (low- or
%   high-pass depending axis and subband, see below). The convolution will
%   save only the central part of the filtered data, i.e. a data length X
%   and kernel length Y will produce a filtered data of length X.
% In sum, the stationary transform does not decimate the obtained
% coefficient-subbands and thus is redundant but protected against
% truncation artifacts.
%
% Inputs:
%
% x             The 3D array to be decomposed.
%
% n             Depth of decomposition.
%
% h             Decomposition low-pass filter.
%
% g             Decomposition high-pass filter.
%
% hz            Optional. If anisotrop WT decomposition is desired then a
%               separate decomposition low-pass filter can be specified for
%               the Z-axis (X & Y will always use filter h and g). If hz is
%               specified then gz must be also specified.
%
% gz            Optional. Decomposition high-pass filter to be used along
%               Z-axis (see hz above).
%
% st            Optional. You can specify whether stationary transform
%               should be used (if st==1) or dyadic (if st==0). If not
%               specified then dyadic WT is used.
%
% outfname      Optional argument. You must precede it with the string
%               argument 'save' to indicate its presence. Name of output
%               .wd3 file. 
%
%
% The output wavelet 3D decomposition structure [c,s] contains the wavelet
% decomposition vector c and the corresponding bookeeping matrix s. 
%   Vector c is organized as:
% c = [ app{n}, h_l{n}, v_l{n}, d_l{n}, a_h{n}, h_h{n}, v_h{n}, d_h{n}, ...
% h_l{n-1}, v_l{n-1}, d_l{n-1}, a_h{n-1}, h_h{n-1}, v_h{n-1}, d_h{n-1}, ...
% h_l{1}, v_l{1}, d_l{1}, a_h{1}, h_h{1}, v_h{1}, d_h{1} ].
%   where the constituents are row vectors such that the low- and high-pass
% filters are applied in different combinations along the X, Y, Z axes. The
% following is the table explaining which filter is applied along the axes
% for each given vector constituent:
%
% NAME X-axis     Y-axis     Z-axis     Comment
% app  low-pass   low-pass   low-pass   approximation coefficients
% h_l  high-pass  low-pass   low-pass   hori. detail coefficients in Z-low
% v_l  low-pass   high-pass  low-pass   vert. detail coefficients in Z-low
% d_l  high-pass  high-pass  low-pass   diag. detail coefficients in Z-low
% a_h  low-pass   high-pass  high-pass  pseudoapproximation coefficients
% h_h  high-pass  low-pass   high-pass  hori. detail coefficients in Z-high
% v_h  low-pass   high-pass  high-pass  vert. detail coefficients in Z-high
% d_h  high-pass  high-pass  high-pass  diag. detail coefficients in Z-high
%
% Each vector is the vectorised storage of a 3D array of coefficients of
% the given subband at a given level (approximation coefficients are
% provided only at the deepest level).
%
% Matrix s is n+2 x 3 matrix such that:
%   s(1,:) = size of app. coef. array (at level n)
%   s(i,:) = size of det. coef. array (at level n-i+2) for i = 2,...,n+1
%   s(n+2,:) = size(x)
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

if ~any((3:9)==nargin)
    error('Wrong number of inputs');
end
if ~any((0:2)==nargout)
    error('Wrong number of outputs');
end
if (ischar(n) || any(n < 1) || any(n ~= fix(n)))
    error('For argument n integer > 0 expected');
end
anisotrop_filter=0;
dwt_st=0;
fid=0;
fpos0=0;
deltaArg=0;
saveWD3=0;
for i=1:length(varargin)
    if strcmp(varargin{i},'save')
        fname=varargin{i+1};
        if isnumeric(fname)
            error('You must specify a filename');
        else
            [~,~,e]=fileparts(fname);
            if strcmp(e,'.wd3')
                saveWD3=1;
                fname2=[fname '1'];
            else
                fname2=fname;
            end
            if saveWD3
                fid=fopen(fname2,'w','l');
                fwrite(fid,1,'double');
                fclose(fid);
            end
            fid=fopen(fname2,'r','l');
            fseek(fid,0,1);
            fpos0=ftell(fid);
            fclose(fid);
        end
        deltaArg=2;
        varargin(i:i+1)=[];
        break;
    end
end
if nargin-deltaArg==3
  [LoF_D,HiF_D] = wfilters(varargin{1},'d'); % you need Wavelet Toolbox for this
end
if nargin-deltaArg==4
  LoF_D = varargin{1};   HiF_D = varargin{2};
end
if nargin-deltaArg==5
  LoF_D = varargin{1};   HiF_D = varargin{2};
  if ischar(varargin{3})
      switch lower(varargin{3})
          case 'st'
              dwt_st=2;
      end
  end
end
if nargin-deltaArg==6
  LoF_D = varargin{1};   HiF_D = varargin{2};
  LoF_D2 = varargin{3};   HiF_D2 = varargin{4};
  anisotrop_filter=1;
end
if nargin-deltaArg==7
  LoF_D = varargin{1};   HiF_D = varargin{2};
  LoF_D2 = varargin{3};   HiF_D2 = varargin{4};
  anisotrop_filter=1;
  if ischar(varargin{5})
      switch lower(varargin{5})
          case 'st'
              dwt_st=2;
      end
  end
end
% Initialization. Precalculate s
s = size(x);
if length(s)~=3
  error('Input must be a 3 dimensional array.');
end
filterlen=repmat(length(LoF_D),1,3);
if anisotrop_filter
    filterlen(3)=length(LoF_D2);
end
for i=1:n
    s=[wapi_wt_trsize(s(1,:), filterlen, dwt_st);s]; 
end
s=s([1 1:end],:);

if saveWD3
    wapi_writevol(fname,s);
    fid0=fopen(fname,'a','l');
    if (fid0<2)
        error('Error opening output file.');
    end
    fwrite(fid0,1,'double');
end
if fid
    [~,lc]=wapi_idxcoef3('c',s,1);
    wapi_extendFile(fname2,lc*8);
    fid=fopen(fname2,'r+','l');
end

c = [];
splitc=0;

if version('-release')<14
    error('Not supported anymore');
else
    dwt3mx=@wapi_mthread_dwt3mex;
    dwt3zmx=@wapi_mthread_dwt3zmex;
end

for i=1:n
    if dwt_st
       if anisotrop_filter
         xd = feval(dwt3zmx,x,16,LoF_D,HiF_D,LoF_D2,HiF_D2,dwt_st); % decomposition
         LoF_D2=dyadup(LoF_D2);
         HiF_D2=dyadup(HiF_D2);
       else
         xd = feval(dwt3mx,x,0,LoF_D,HiF_D,dwt_st); % decomposition
       end

       LoF_D=dyadup(LoF_D);
       HiF_D=dyadup(HiF_D);

    else
       if anisotrop_filter
         xd = feval(dwt3zmx,x,0,LoF_D,HiF_D,LoF_D2,HiF_D2); % decomposition
       else
         xd = feval(dwt3mx,x,0,LoF_D,HiF_D); % decomposition
       end
    end
    if fid
        f=wapi_idxcoef3('c',s,i);
        fseek(fid,fpos0+(f-1)*8,-1);
        c(i)=0; %#ok<*AGROW>
        for b=2:8
            c(i)=c(i)+fwrite(fid,xd{b},'double');
        end
        xd(2:8)=[];
        x = xd{1};
        clear xd
    else
        x = xd{1};
        xd{1}=[];
        if splitc
            c{i}=xd; %#ok<UNRCH>
        else
            for j=8:-1:2
                c=[xd{j}(:)' c]; % store details
            end
        end
    end
    sz=size(x);
    if ~isequal(sz,s(end-i,:))
        error('Level %d details have incorrect size ([%s] vs. [%s])',num2str(sz),num2str(s(end-i,:)));
    end
%    s = [size(x);s];                 %#ok<AGROW> % store size
    fprintf(1,'Level %d done.\n',i);
end

% Last approximation.
clear xd LoF_D HiF_D LoF_D2 HiF_D2
if fid
    fseek(fid,fpos0+(wapi_idxcoef3('app',s,1)-1)*8,-1);
    c(end+1)=fwrite(fid,x(:),'double');
    fclose(fid);
    if saveWD3
        fwrite(fid0,sum(c),'double');
        fclose(fid0);
    end        
else
    if splitc
        c{i+1}=x; %#ok<UNRCH>
    else
        c = [x(:)' c];
    end
end
sz=size(x);
if ~isequal(sz,s(1,:))
    error('Approximation has incorrect size ([%s] vs. [%s])',num2str(sz),num2str(s(1,:)));
end

function y=dyadup(x)

l = 2*length(x)+1;
y = zeros(1,l);
y(2:2:l) = x;
