function voldat=wapi_getframe(petvol, fr)
%wapi_getframe          Read 3D image data for a single frame
%
% voldat=wapi_getframe(petvol, fr)
%
% Read the image data of a single frame into memory. The image header
% information must already have been read into memory.
%
% Inputs:
%
% petvol        A MATLAB struct describing a PET image as provided by
%               <a href="matlab: help wapi_readim">wapi_readim</a>.
%
% fr            The frame number to be read (starts at 1).
%
%
% Outputs:
%
% voldat        3D array of quantified voxel data for the frame.
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

switch petvol.type
    case 'nii'
        voldat=petvol.private.dat(:,:,:,fr).*petvol.frameinfo(fr,3);
	case 'nis'
		voldat=petvol.private(fr).dat(:,:,:).*petvol.frameinfo(fr,3);
end
