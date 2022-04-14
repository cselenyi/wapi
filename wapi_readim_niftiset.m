function imvol=wapi_readim_niftiset( nisname )
%wapi_readim_niftiset      Read image info for set of 3D NIfTI-1 images
%
% imvol=wapi_readim_niftiset( niiname )
%
% Read image header (volume & timing) information in multi-file
% NIfTI-1 format. SPM5 or later must be installed in order to be able to
% read set of NIfTI-1 images.
%
% Inputs:
%
% infname       Virtual name of the input NIfTI-1 file set. The actual 3D
%               frame filenames should be named (after dropping the ".nis"
%               extension from infname) infname_001.nii, infname_002.nii,
%               ... infname_NNN.nii, i.e. frames should be numbered with
%               leading zeros so that they line up correctly when sorted
%               alphanumerically! The 3D frame images should contain timing
%               information: the start of the frame should be in the
%               ".timing.toffset" field and the frame duration in the
%               ".timing.tspace" field (both in seconds) of the NIFTI
%               object obtained when opening the frame image using "nifti".
%               The frames should be correctly scaled to radioactivity
%               units.
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

imvol.filename=nisname;

imvol.type='nis';

[pth,basename]=fileparts(nisname);
d=dir(fullfile(pth,[basename '_*.nii']));
frameNames=sort({d.name});
numFr=length(frameNames);
removeId=false(numFr,1);
fr=1;
for i=1:numFr
	frName=frameNames{i};
	toks=regexp(frName,'(.+)_(\d+)\.nii','tokens');
	if isempty(toks) || ~strcmp(toks{1}{1},basename)
		removeId(i,1)=true; % this is some other file, let's remove it in frameinfo
		continue;
	end
	frNum=str2double(toks{1}{2});
	if frNum~=fr
		error('Frame %d has a filename %s with a non-matching frame number %d', fr, frName, frNum);
	end
	if fr==1
		N=nifti(fullfile(pth,frName));
		imvol.size=N.dat.dim;
		if length(imvol.size)==4 && imvol.size(4)>1
			error('Frame %d is a 4D image, not a 3D one', fr);
		end
		imvol.mat=N.mat;
		imvol.frameinfo=zeros(numFr,3);
	else
		N(fr)=nifti(fullfile(pth,frName));
		if ~isequal(imvol.size(1:3),N(fr).dat.dim(1:3))
			error('Frame %d array dimensions do not match that of frame one', fr);
		end
		if ~isequal(imvol.mat,N(fr).mat)
			error('Frame %d space definition does not match that of frame one', fr);
		end
	end
	try
		imvol.frameinfo(i,:)=[N(fr).timing.toffset N(fr).timing.tspace 1]; % 3 columns: frame start, frame duration in sec and frame scaling factor
	catch
		error('Could not read frame timing information for frame %d', fr);
	end
	fr=fr+1;
end
imvol.frameinfo(removeId,:)=[];
imvol.size(4)=size(imvol.frameinfo,1);
imvol.private=N; % we cannot directly place it in private, it spits an error