function wapi_dwt3dyn(infname,outfname,d,varargin)
%wapi_dwt3dyn           Wavelet transform 4D PET image frame-by-frame
%
% wapi_dwt3dyn(infname,outfname,d,h,g)
% wapi_dwt3dyn(infname,outfname,d,h,g,hz,gz)
% wapi_dwt3dyn(infname,outfname,d,h,g,st)
% wapi_dwt3dyn(infname,outfname,d,h,g,hz,gz,st)
% wapi_dwt3dyn(infname,outfname,d,...,'limits',lims,...)
%
% Transforms the input dynamic (4D) PET volume to wavelet space one frame
% at a time using the given decomposition low- and highpass filters and
% specified maximum depth of decomposition. The wavelet transform (WT) can
% be limited in space as time (see help below).
% Saves output into a '.wd3' file (3-d dynamic wavelet transform file).
%
% Inputs:
%
% infname       Name of input 4D PET image.
%
% outfname      Name of output .wd3 file. Each frame is saved to a separate
%               file in case of stationary transform (e.g. to first frame
%               to outfname.wd31, 2nd to *.wd32, etc.).
%
% d             The maximum depth of decomposition.
%
% h             Decomposition low-pass filter. See <a href="matlab: help uvi_lemarie">uvi_lemarie</a>.
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
% lims          Optional argument. You must precede it with the string
%               argument 'limits' to indicate its presence. This argument
%               can be used to limit the WT to only certain frames and/or a
%               certain sub-volume of the PET image. See help for
%               <a href="matlab: help wapi_wapi">wapi_wapi</a> for a detailed description.
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
%     WAPI 1.1 2022-04-14


nArgs=nargin;
limits=[];
for i=1:numel(varargin)
    if ischar(varargin{i}) && strcmpi(varargin{i},'limits')
        limits=varargin{i+1};
        varargin(i:i+1)=[];
        nArgs=nArgs-2;
        break;
    end
end
if nArgs==5 || nArgs==7
  st=0;
end
if nArgs==6 || nArgs==8
  if nArgs==6
    st=varargin{3};
  else
    st=varargin{5};
  end
end
if exist('st') ~= 1 %#ok<EXIST>
  error('Incorrect number of inputs.');
end
H=varargin{1};
G=varargin{2};
anisotrop_filter=0;
if nArgs>6
  H2=varargin{3};
  G2=varargin{4};
  anisotrop_filter=1;
end
if ischar(infname)
    invol=wapi_readim(infname, tmsname);
else
    if isstruct(infname)
        invol=infname;
        infname=invol.filename;
    else
        error('1st input must be filename or PET volume structure from wapi_readim_...');
    end
end
if isempty(invol)
  error('Unable to load PET file ''%s''',infname);
end
if invol.size(4)==1
  error('Input volume is not dynamic.');
end
sz=invol.size;

if ~isempty(limits)
    framelimits=limits(:,4)';
else
    framelimits=[1 sz(4)];
end
for fr=framelimits(1):framelimits(2)
    volume_dat=wapi_getframe(invol,fr);
    if ~isempty(limits)
        volume_dat=volume_dat(limits(1,1):limits(2,1),limits(1,2):limits(2,2),limits(1,3):limits(2,3));
    end
    fr_fname=strcat(outfname,num2str(fr-framelimits(1)+1));
    fid=fopen(fr_fname,'w','l');
    fwrite(fid,fr-framelimits(1)+1,'double');
    fclose(fid);
    if st
        if anisotrop_filter
            [c,s]=wapi_wavedec3st(volume_dat,d,H,G,H2,G2,'st','save',fr_fname);
        else
            [c,s]=wapi_wavedec3st(volume_dat,d,H,G,'st','save',fr_fname);
        end
    else
        if anisotrop_filter
            [c,s]=wapi_wavedec3st(volume_dat,d,H,G,H2,G2,'save',fr_fname);
        else
            [c,s]=wapi_wavedec3st(volume_dat,d,H,G,'save',fr_fname);
        end
    end
    drawnow
    if fr==framelimits(1)
        wapi_writevol(outfname,s);
        fid=fopen(outfname,'a','l');
        if (fid<2)
            error('Error opening output file.');
        end
        fwrite(fid,diff(framelimits)+1,'double');
        fwrite(fid,sum(c),'double');
        fclose(fid);
    end
    fprintf(1,'Frame %d decomposed and saved.\n',fr);
    drawnow
    clear c s
end
return;
