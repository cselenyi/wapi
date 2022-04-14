function out=wapi_readim(infname, tmsname)
%wapi_readim            Read 3D/4D (PET) image info.
%
% imvol=wapi_readim(infname, [tmsname])
%
% Read 3D/4D (PET) image header in one of the supported file formats:
% NIfTI-1 (.nii)
%
% Inputs:
%
% infname       Name of the input (4D) PET image file.
%
% tmsname       Name of extra file needed for certain formats which
%               contains frame timing and/or scaling information. Currently
%               needed for:
%               Single-file 4D NIfTI format: see help for
%               <a href="matlab: help wapi_readim_nifti">wapi_readim_nifti</a> for exact description of the extra file.
%
%
% Outputs:
%
% imvol         A MATLAB structure containing the header data for the PET
%               image in a format required by wapi_wapi.
%
%
% Invoking wapi_readim without input arguments will list the available PET
% format readers in a MATLAB structure.
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

persistent formats

if isempty(formats)
    formats.nii=struct( ...
        'name', 'nifti', ...
        'desc', '4-D NIfTI-1 image', ...
        'tms', 1);
    formats.nis=struct( ...
        'name', 'niftiset', ...
        'desc', 'Set of 3-D NIfTI-1 images (1st frame)', ...
        'tms', 0);
%     formats.vff=struct( ...
%         'name', 'vff', ...
%         'desc', 'VFF images', ...
%         'tms', 0);
end

if ~nargin
    out=formats;
    return;
end

exts=fieldnames(formats);
[inPath,inName,inExt]=fileparts(infname);
inExt(1)=[];
if ~any(strcmpi(exts, inExt))
    error('Cannot handle PET images of type .%s', inExt);
end

format=formats.(inExt);
if format.tms && (nargin<2 || isempty(tmsname))
    error('PET images of type .%s require additional input file', inExt);
end

funName=sprintf('wapi_readim_%s', format.name);

if format.tms
    invol=feval(funName,infname, tmsname);
else
    invol=feval(funName,infname);
end

out=invol;