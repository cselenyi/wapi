function imvol=wapi_readim_nifti( niiname, tmsname )
%wapi_readim_nifti      Read 3D/4D NIfTI-1 image info
%
% imvol=wapi_readim_nifti( niiname, tmsname )
%
% Read (4D) PET image header (volume & timing) information in single-file
% NIfTI-1 format. SPM5 or later must be installed in order to be able to
% read NIfTI-1 images.
%
% Inputs:
%
% infname       Name of the input (4D) PET NIfTI-1 file.
%
% tmsname       Name of extra file which contains frame timing (and
%               scaling) information for 4D PET image data. The timing file
%               should contain no header lines and:
%               A single column of values indicating frame middle times in
%               seconds
%                 OR
%               Two columns of values indicating frame starting and
%               duration times in seconds.
%                 OR
%               Three columns of values indicating frame starting and
%               duration times in seconds and frame scaling factors. When
%               image data is read into memory the scaling factor is
%               applied to obtain the final quantified voxel values.
%
%               If omitted then fake image timing information is created.
%
% Outputs:
%
% imvol         A MATLAB structure containing the header data for the PET
%               image in a format required by wapi_wapi.
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

% imvol.filename
% imvol.type
% imvol.size
% imvol.mat
% imvol.frameinfo
% imvol.private

imvol.filename=niiname;

imvol.type='nii';

N=nifti(niiname);

imvol.size=N.dat.dim;
if length(imvol.size)==3
    imvol.size(4)=1;
end
imvol.mat=N.mat;

if nargin>1 && ~isempty(tmsname)
    tms=textread(tmsname);
    switch size(tms,2)
        case 1 % only frame centers provided
            dr=diff([0;tms]);
            st=tms-dr/2;
            imvol.frameinfo=[st dr ones(imvol.size(4),1)]; % 3 columns: frame start, frame duration in sec and frame scaling factor
        case 2 % frame start and duration provided
            imvol.frameinfo=[tms ones(imvol.size(4),1)];
    end
else
    imvol.frameinfo=[(0:imvol.size(4)-1)' repmat([eps 1],imvol.size(4),1) ones(imvol.size(4),1)];
end
imvol.private=N;
