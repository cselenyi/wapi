function outfnames=wapi_wapi(infname,tmsname,reffname,ncoeff,depth,stationary,outfname,numpoints,options) % ,refmaskfile,weights,k2ref,lims)
%wapi_wapi              Wavelet Aided Parametric Imaging (WAPI)
%
% Wavelet estimation of dynamic (4D) PET images using 3D wavelet analysis
% and Logan's plot (using either reference region or plasma input
% function). The decomposition is done using Battle-Lemarie wavelet
% filters.
%   The input 4D PET image is assumed to be decay-corrected. Multi-linear
% (reference-)Logan fit is used to calculate a parametric image of Total
% distribution volume (V_T) or Distibution Volume Ratio (DVR) depending on
% the input function used.
%   The wavelet transform (WT) data of the original PET image (and of the
% reference mask if it is used) is deleted after the WAPI processing is
% finished (these are created in the temporary directory). However, the WT
% of the parametric image itself (as well as the WT of the standard
% deviation of the Logan fit, which is not inverted back to image space) is
% retained (this file is created in the same folder as the output
% parametric image).
%
% Version 2022-09-28
%
% Usage:
% wapi_wapi(infname,tmsname,reffname,ncoeff,depth,stationary,outfname, ...
%              numpoints,[options])
%
% Inputs:
%
%    infname        <a href="matlab: try;wapi_gui;uicontrol(findobj('Tag','infnameEdit'));catch;end;">4D PET image</a>:
%                   Name of the file containing the input 4D PET image
%                   data. See help for <a href="matlab: help wapi_readim">wapi_readim</a> for a list of supported
%                   file types.
%
%    tmsname        <a href="matlab: try;wapi_gui;uicontrol(findobj('Tag','tmsnameEdit'));catch;end;">Frame (scaling &) timing file</a>:
%                   Name of additional file containing PET frame timing
%                   and/or scaling information. Required only in case of
%                   certain 4D PET image file types.
%
%    reffname       <a href="matlab: try;wapi_gui;uicontrol(findobj('Tag','reffnameEdit'));catch;end;">Input function TAC file</a>:
%                   Reference region radioactivity file or blood (plasma)
%                   input function file. See help for <a href="matlab: help wapi_readref">wapi_readref</a> for
%                   expected format of reference-region input data. See
%                   help for <a href="matlab: help wapi_readblood">wapi_readblood</a> for expected format of plasma
%                   input data (the timing in the blood file is assumed to
%                   be in seconds). In order to indicate that the input is
%                   from the blood (plasma) the name of the file either has
%                   to contain the word 'blood' or the extension should be
%                   '.blo'.
%
%    ncoeff         <a href="matlab: try;wapi_gui;uicontrol(findobj('Tag','ncoeffEdit'));catch;end;">Length of wavelet filter kernel (# of coeffs)</a>:
%                   Number of coefficients of Battle-Lemarie wavelet
%                   filters (it should be a positive even number in the
%                   range of 4 and 128). Typical value for HRRT PET data is
%                   16.
%
%    depth          <a href="matlab: try;wapi_gui;uicontrol(findobj('Tag','depthEdit'));catch;end;">Depth of wavelet decomposition</a>:
%                   Number of decomposition levels to calculate. Value
%                   should be in the range of 1 and 8 although above 5 it
%                   appears impractical. Recommended value for HRRT PET
%                   data is 3.
%
%    stationary     <a href="matlab: try;wapi_gui;uicontrol(findobj('Tag','stationaryCheckbox'));catch;end;">Stationary (non-dyadic) wavelet transform</a>:
%                   Flag indicating whether stationary wavelet transform
%                   (true) or dyadic stransform (false) should be used. The
%                   stationary wavelet transform (WT) is redundant so at
%                   each iteration level of the WT 7 times the size of the
%                   original PET data (in double-precision floating point
%                   representation) is added to the decomposition
%                   structure for the details coefficients and at the final
%                   iteration level the approximation coeffcients require
%                   also the same size as the original PET frame data. So
%                   for example with a WT of depth 3 the whole
%                   decomposition structure contains 22 times (7*3+1) the
%                   size of the original PET frame data (e.g. in case of
%                   HRRT data it is 2.22 GB of storage for each frame).
%                     If using the dyadic WT the coefficients are decimated
%                   so the transform is non-redundant and thus it will
%                   always equal the size of the original PET image data
%                   (albeit in double-precision floating point format). In
%                   case of HRRT it means that only 103.5 MB is required
%                   for each frame at any decomposition depth. The drawback
%                   of the dyadic WT is that it can lead to artifacts
%                   especially in case of noisier data or high-contrast
%                   image features (e.g. in case of radioligands with high
%                   specific binding in striatum). Therefore the stationary
%                   transform is superior and is recommended if
%                   computational resources allow. Tip: use a fast SSD or
%                   ramdisk for your temporary storage. Set the TMP
%                   environment variable in Linux or the TEMP environment
%                   variable in Windows to the fast storage location when
%                   starting up MATLAB and before using WAPI. There can be
%                   tremendous gain in processing speed!
%
%    numpoints      <a href="matlab: try;wapi_gui;uicontrol(findobj('Tag','numpointsEdit'));catch;end;">Number of points fitted on (reference-)Logan plot</a>:
%                   Number of last points on the (ref.)Logan plot that are
%                   fitted to obtain the kinetic parameter (V_T or DVR).
%
%    options        An optional matlab structure specifying additional
%                   parameters. The following fields (options) are
%                   recognized:
%
%       refmaskfile <a href="matlab: try;wapi_gui;uicontrol(findobj('Tag','refmaskfileEdit'));catch;end;">3D mask of reference region</a>:
%                   Optional. 3D-mask with reference region (filename of
%                   .mat or .nii file or variable in memory with 3D data).
%                   If it is a .mat file then the (alphabetically) first
%                   variable in the .mat file should be the 3D mask data.
%                   Default is empty: [].
%                     The 3D mask is also transformed into wavelet space
%                   and for each details subband this mask-transform is
%                   used to obtain a threshold value in the following way:
%                     1. The mean of the absolute value of all those
%                     coefficients in the given subband of the
%                     mask-transform is calculated that are not equal to
%                     zero. 
%                     2. The index of all coefficients in the given subband
%                     of the mask-transform is stored that have an absolute
%                     value greater than the mean value calculated in the
%                     previous step. These are the "reference-region"
%                     coefficients.
%                     3. The area-under-the-curve (AUC) of each
%                     coefficient-time curve of the given subband of the WT
%                     of the PET image data is calculated.
%                     4. The mean absolute AUC is calculated for the 
%                     reference-region coefficients.
%                     5. The threshold value is 1% of the mean value. Any
%                     coefficient that has a lower absolute AUC than the
%                     threshold is excluded from the Logan fit and the
%                     output parameter is set to 0.
%                     If the reference mask is not provided (which is
%                   always the case if plasma input is used) then no
%                   thresholding is performed (calculation is done for all
%                   coefficients).
%
%       weights     <a href="matlab: try;wapi_gui;uicontrol(findobj('Tag','weightsEdit'));catch;end;">Weights</a>:
%                   Optional. Frame-by-frame weighting factors to be used
%                   in the kinetic analysis. Name of a text file containing
%                   one number per line of factors or a vector. If omitted
%                   then equal weights are assumed.
%
%       k2ref       <a href="matlab: try;wapi_gui;uicontrol(findobj('Tag','k2refEdit'));catch;end;">k2 in reference region</a>:
%                   Optional. The k2 rate constant in the reference region.
%                   If omitted then 0 is assumed (i.e. no correction). Of
%                   course it should be 0 (i.e. omitted) in case plasma
%                   input.
%
%       lims        <a href="matlab: try;wapi_gui;uicontrol(findobj('Tag','limsEdit'));catch;end;">Spatiotemporal limits</a>:
%                   Optional. Limit the WT and thus the calculations to
%                   part of the image (in  space and/or time).
%                   It can be a matrix with size:
%                     2 x 3 for [xmin, ymin, zmin; xmax, ymax, zmax]
%                   OR
%                     2 x 4 for [xmin, ymin, zmin frmin; xmax, ymax, zmax, frmax]
%                   OR
%                     vector of length 2 for [frmin frmax]
%                   Limits are inclusive and spatial limits are defined in
%                   voxel coordinates starting at 1 (e.g. in case of HRRT
%                   they should be 1-256 for x and y and 1-207 for z axes,
%                   respectively).
%
%        targetmaskfile <a href="matlab: try;wapi_gui;uicontrol(findobj('Tag','targetmaskfileEdit'));catch;end;">3D mask of target region</a>:
%                   Optional. Limit the WT to only certain region of the
%                   image. Voxels outside the target mask will be zeroed
%                   before the wavelet transform. Useful for suppressing
%                   artifacts in the parametric image from high
%                   radioactivity voxels outside the region of interest,
%                   for example whole brain. Should be the name of a mask
%                   file specifying voxels to be included in the analysis.
%
%        targetmaskpadding <a href="matlab: try;wapi_gui;uicontrol(findobj('Tag','targetmaskpaddingEdit'));catch;end;">Padding of target mask</a>:
%                   Optional. The target mask can be padded, i.e. dilated
%                   to increase its coverage at its borders. It should be
%                   the number of voxels to expand the border with.
%
%        deleteoutputwd3 <a href="matlab: try;wapi_gui;uicontrol(findobj('Tag','deleteoutputwd3Checkbox'));catch;end;">Delete output WT file</a>:
%                   Optional. If set (checked) then the WT file of the
%                   calculated parametric image will be deleted (but not
%                   the actual image file or course :-) ).
%
%  Outputs:
%
%    outfname       <a href="matlab: try;wapi_gui;uicontrol(findobj('Tag','outfnameEdit'));catch;end;">Output 3D parametric image</a>:
%                   Parametric image is saved in .nii file with this base
%                   name. The number of coeffs and depth of decomposition,
%                   and the type of fitting used ('mllog' for multi-linear
%                   Logan) is inserted into the filename. 
%                     E.g.: dymmy_wapi.nii will become:
%                   dummy_wapi16-3_mllog.nii with a wavelet transform with
%                   kernel length of 16 and decomposition depth of 3, and 
%                   'mllog' parameter estimator (multi-linear Logan fit).
%
%  Example:
%
%  wapi_wapi('dummy.nii','dummy_times.txt','dummy_refTAC.txt',16,3,1,
%            'dummy_wapist.nii',6,[],'dummy_weights.txt',0.15)
%
%      Output will be: dummy_wapist16-3_mllog.nii with its wavelet transform
%      in 'dummy_wapist16-3_mllog_p1.wd3*' and the WT of the standard
%      deviation of the fit in 'dummy_wapist16-3_mllog_std.wd3*'. 
%  
% Please use the following reference in publications based on WAPI results:
%
%<a href="http://www.ncbi.nlm.nih.gov/pubmed/16859930">Cselenyi Z, Olsson H, Halldin C, Gulyas B, Farde L.:
%A comparison of recent parametric neuroreceptor mapping approaches based
%on measurements with the high affinity PET radioligands [11C]FLB 457 and
%[11C]WAY 100635. Neuroimage. 2006 Oct 1;32(4):1690-708.</a>
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
%     WAPI 1.2 2022-09-27

if nargin<8 || nargin>9
  error('Number of inputs must be 8 - 9. See help.');
end

defOptions=struct('refmaskfile',[],'weights',[],'k2ref',[],'lims',[],'targetmaskfile','','targetmaskpadding',1,'deleteoutputwd3',false);
if nargin>6
	if ~isstruct(options)
		error('7th input must be a struct. See help.');
	end
	opts=options;
	options=defOptions;
	fns=fieldnames(options);
	for f=1:length(fns)
		fld=fns{f};
		if isfield(opts,fld)
			options.(fld)=opts.(fld);
		end
	end
else
	options=defOptions;
end
fns=fieldnames(options);
for f=1:length(fns)
	fld=fns{f};
	eval(sprintf('%s=options.%s;',fld,fld));
end

framelimits=[];
if ischar(numpoints)
    numpoints=str2num(numpoints); %#ok<ST2NM>
    if isempty(numpoints)
        error('numpoints must be a scalar or vector!');
    end
end

if isempty(lims) %#ok<NODEF>
    lims=framelimits;
else
    if ischar(lims)
        lims=str2num(lims); %#ok<ST2NM>
        if isempty(lims)
            error('lims must be a numeric array');
        end
    end
end
if isempty(k2ref) %#ok<NODEF>
    k2ref=0.0; % 0.1;
elseif ischar(k2ref)
    k2ref=str2double(k2ref);
    if isnan(k2ref)
        error('k2ref must be a number if given');
    end
end
if ischar(weights) && ~isempty(weights) %#ok<NODEF>
	weights=textread(weights); %#ok<DTXTRD> % it should be a simple file with one weighting value per line
end
plasmainput=0;
if ischar(reffname)
    if ~isempty(strfind(reffname,'.blo')) || ~isempty(strfind(reffname,'blood'))
        ref=wapi_readblood(reffname,'min');
        ref=ref(:,[1 3]); % throw away total blood
        plasmainput=1;
    else
        ref=wapi_readref(reffname);
    end
    if plasmainput
        disp('Plasma input function detected. Output will be V_T.');
    end
else
    ref=reffname;
end

if ischar(ncoeff)
    ncoeff=str2double(ncoeff);
    if isnan(ncoeff)
        error('ncoeff must be a number!');
    end
end
[h,g,rh,rg]=uvi_lemarie(ncoeff); % Battle-Lemarie wavelet decomposition and reconstruction filters 
ncoeff=length(h);
if ischar(depth)
    depth=str2double(depth);
    if isnan(depth)
        error('depth must be a number!');
    end
end

hz=[]; %#ok<NASGU>
if ~ischar(infname)
    error('1st input must be a valid PET filename');
end
[inPath,inName]=fileparts(infname); %#ok<ASGLU>
invol=wapi_readim(infname, tmsname);

invol.orisize=invol.size(1:3);
sz=invol.size;
time=(invol.frameinfo(:,1)./60) + (invol.frameinfo(:,2)./60./2);
if ~plasmainput
    ref=[time ref];
end
if isempty(weights)
    fprintf(1,'Uniform weights assumed.\n');
    drawnow
    weights=ones(sz(4),1);
end

Z = spm_imatrix(invol.mat);
asp = Z(7:9); % voxel size
R = spm_matrix([0 0 0 Z(4:6)]); % rotational components
R = R(1:3,1:3); % direction cosines
[v,ipermOrder]=max(abs(R)); %#ok<ASGLU> % where x,y,z dim appears in input spm image
if ~isequal(ipermOrder, [1 2 3])
    error('Input PET data must be stored in x-y-z arrays');
end

if ~isempty(targetmaskfile) && ischar(targetmaskfile)
	% a file, load it
	[~,~,ext]=fileparts(targetmaskfile);
	switch ext
		case '.mat' % target mask is in a MATLAB .mat file (should be 1st variable in .mat as a 3D array)
			str=load(targetmaskfile);
			nm=fieldnames(str);
			targetmask=str.(nm{1});
			if ~isequal(size(targetmask),invol.orisize)
				error('First variable of target mask .mat file must be of the same size as input image.');
			end
		case '.nii'
			targetmask=wapi_readim_nifti(targetmaskfile);
			targetmask_mat=targetmask.mat;
			targetmask=wapi_getframe(targetmask,1);
			if ~isequal(size(targetmask),invol.orisize)
				error('Target mask must be of the same array size as input image.');
			end
			im2mMat=inv(invol.mat)*targetmask_mat; %#ok<MINV>
			cosdiff=abs(abs(im2mMat(1:3,1:3))-eye(3));
			if any(cosdiff(:)>1e-10)
				error('Voxel layout of target mask must be the same as that of the input image, except for flipping.');
			end
			for d=1:3
				if im2mMat(d,d)<0
					if im2mMat(d,4)~=invol.orisize(d)+1
						error('Voxel layout of target mask must be the same as that of the input image, except for flipping.');
					end
					targetmask=flip(targetmask,d);
				elseif im2mMat(d,4)~=0
					error('Voxel layout of target mask must be the same as that of the input image, except for flipping.');
				end
			end
			targetmask(targetmask<.01)=0;
			targetmask(targetmask>=.01)=1;
		otherwise
			error('Unsupported mask filetype %s',ext);
	end
elseif ~isempty(targetmaskfile) && (isnumeric(targetmaskfile) || islogical(targetmaskfile)) % just an array
	targetmask=targetmaskfile;
	if ~isequal(size(targetmask),invff.orisize.dim(1:3))
		error('Target mask must be of the same size as input image.');
	end
else
	targetmask=[];
end
if ~isempty(targetmask)
	targetmask=targetmask>0;
	if targetmaskpadding>0
		targetmask=imdilate(targetmask,strel_bol(targetmaskpadding));        % erode voxels at boundary in MR space
	end
	fprintf(1,'Will use target mask for calculating wavelet transform.\n');
	drawnow
end

if ~isempty(lims)
    switch size(lims,2)
        case 4
            if size(lims,1) ~=2
                error('Limits specified erroneously!');
            end
        case 3
            if size(lims,1) ~=2
                error('Limits specified erroneously!');
            end
            lims(:,4)=[1; sz(4)];
        case 2
            if size(lims,1) ~=1
                error('Limits specified erroneously!');
            end
            lims=[1 1 1 lims(1);sz(1) sz(2) sz(3) lims(2)];
        case 1
            if size(lims,1) ~= 2
                error('Limits specified erroneously!');
            end
            lims=[1 1 1 lims(1);sz(1) sz(2) sz(3) lims(2)];
        otherwise
            error('Limits specified erroneously!');
    end
    sz=lims(2,:)-lims(1,:)+1;
    invol.size=sz; %(1:3);
    ref=ref(ref(:,1)>=time(lims(1,4))-10*eps & ref(:,1)<=time(lims(2,4))+10*eps,:);
    time=time(lims(1,4):lims(2,4));
    weights=weights(lims(1,4):lims(2,4));
	if ~isempty(targetmask)
		targetmask=targetmask(lims(1,1):lims(2,1),lims(1,2):lims(2,2),lims(1,3):lims(2,3),:);
	end
end

factor=asp(1) / asp(3);
l2=round(ncoeff * factor);
if mod(l2,2)==1
    l2=l2+1; %#ok<NASGU>
end
l2=ncoeff;
fprintf(1,'Kernel lengths: xy: %d   z: %d\n',ncoeff, l2);
drawnow
[hz, gz, rhz, rgz]=uvi_lemarie(l2);
wd3name=fullfile(tempdir,[inName '.wd3']);
doit=1;
if exist(wd3name,'file')==2 
    fprintf(1,'Wavelet transform file ''%s'' already exists.\n',wd3name);
    drawnow
    while 1
        r='y'; %input('Do you want to use it [Y/N]? ', 's');
        if ~( strcmpi( r, 'y' ) || strcmpi( r, 'n') )
            beep;
            continue;
        end
        if strcmpi( r, 'y' )
            doit=0;
        end
        break;
    end
end
if doit
    fprintf(1,'Converting ''%s'' from normal to wavelet space in ''%s''\n',infname,wd3name);
    drawnow
end
if stationary
    fprintf(1,'Stationary wavelet transform is used.\n');
    drawnow
    if doit
        wapi_dwt3dyn(invol,wd3name,depth,h,g,hz,gz,1,'limits',lims,'targetmask',targetmask);
    end
else
    if doit
        wapi_dwt3dyn(invol,wd3name,depth,h,g,hz,gz,0,'limits',lims,'targetmask',targetmask);
    end
end
if doit
    fprintf(1,'ready.\n');
    drawnow
end

if ischar(outfname)
  odp=strfind(outfname,'.nii');
  if isempty(odp)
    error('Outfilename must have extension ''.nii''.');
  end
else
  error('7th input must be a filename.');
end
 
if any(size(ref(:,1))~=size(time)) % resample input function to PET timeframes but keep peak
    ref=wapi_matchIF2Frames(ref, time);
    if any(size(ref(:,1))~=size(time)) || any(ref(:,1)~=time)
        interpolate=1;
        petTime=time;
        time=ref(:,1);
        weights=interp1(petTime,weights,time);
        disp('Voxel/coefficient TACs will be interpolated to frame-matched input function!');
    else
        interpolate=0;
        petTime=time;
    end
else
    interpolate=0;
    petTime=time;
    ref(:,1)=time;
end

[offs,numfr,lc,s]=wapi_readwd3(wd3name); %#ok<ASGLU>
if ~isempty(lims)
    numfr=sz(4);
end

userefmask = 1; % set to zero to quickly disable this part of the code
if userefmask
    if ~isempty(refmaskfile) && ischar(refmaskfile)
        [~,~,ext]=fileparts(refmaskfile);
        switch ext
            case '.mat' % mask is in a MATLAB .mat file (should be 1st variable in .mat as a 3D array)
                wd3mask=fullfile(tempdir,[strrep(strrep(strrep(strrep(refmaskfile,'.mat',''),'\','_'),'/','_'),':','') '.wd3']);
                str=load(refmaskfile);
                nm=fieldnames(str);
                refmask.volume=str.(nm{1});
                if ~isequal(size(refmask.volume),invol.orisize)
                    error('First variable of mask .mat file must be of the same size as input image.');
                end
                if ~isempty(lims)
                    refmask.volume=refmask.volume(lims(1,1):lims(2,1),lims(1,2):lims(2,2),lims(1,3):lims(2,3),:);
                end
            case '.nii'
                refmask=wapi_readim_nifti(refmaskfile);
				wd3mask=fullfile(tempdir,[strrep(strrep(strrep(strrep(refmaskfile,'.nii',''),'\','_'),'/','_'),':','') '.wd3']);
                refmask.volume=wapi_getframe(refmask,1);
				if ~isequal(size(refmask.volume),invol.orisize)
					error('Reference mask must be of the same array size as input image.');
				end
				im2mMat=inv(invol.mat)*refmask.mat; %#ok<MINV>
				cosdiff=abs(abs(im2mMat(1:3,1:3))-eye(3));
				if any(cosdiff(:)>1e-10)
					error('Voxel layout of reference mask must be the same as that of the input image, except for flipping.');
				end
				for d=1:3
					if im2mMat(d,d)<0
						if im2mMat(d,4)~=invol.orisize(d)+1
							error('Voxel layout of reference mask must be the same as that of the input image, except for flipping.');
						end
						refmask.volume=flip(refmask.volume,d);
					elseif im2mMat(d,4)~=0
						error('Voxel layout of reference mask must be the same as that of the input image, except for flipping.');
					end
				end
				refmask.volume(refmask.volume<.01)=0;
                refmask.volume(refmask.volume>=.01)=1;
                if ~isempty(lims)
                    refmask.volume=refmask.volume(lims(1,1):lims(2,1),lims(1,2):lims(2,2),lims(1,3):lims(2,3),:);
                end
            otherwise
                error('Unsupported mask filetype %s',ext);
        end
        mh=h;
        mg=g;
        disp('Creating reference mask in wavelet space.');
        if stationary
            [cm,sm]=wapi_wavedec3st(refmask.volume,depth,mh,mg,hz,gz,'st','save',wd3mask);  %#ok<ASGLU>
        else
            [cm,sm]=wapi_wavedec3st(refmask.volume,depth,mh,mg,hz,gz,'save',wd3mask);  %#ok<ASGLU>
        end
        clear sm refmask
    else
        if ~isempty(refmaskfile) && isqual(size(refmaskfile),invol.orisize) % mask provided as a 3D array
            wd3mask=[tempname '_mask.wd3'];
            mh=h;
            mg=g;
            refmask=refmaskfile;
            if ~isempty(lims)
                refmask=refmask(lims(1,1):lims(2,1),lims(1,2):lims(2,2),lims(1,3):lims(2,3),:);
            end
            disp('Creating reference mask in wavelet space.');
            if stationary
                [cm,sm]=wapi_wavedec3st(refmask,depth,mh,mg,hz,gz,'st','save',wd3mask);  %#ok<ASGLU>
            else
                [cm,sm]=wapi_wavedec3st(refmask,depth,mh,mg,hz,gz,'save',wd3mask);  %#ok<ASGLU>
            end
            clear sm
        else
            cm=[];%#ok<NASGU> 
            userefmask=0;
        end
    end
end % if userefmask==1

% Modelling
fprintf(1,'Kinetic Analysis:\n');
drawnow

timeL=time.*60; % seconds
inputL=[ref(:,1)*60 ref(:,2)];

plus=1; %#ok<NASGU>

cnames={'h_l' , 'v_l' , 'd_l' , 'a_h' , 'h_h' , 'v_h' , 'd_h'  };
for lev=1:depth
    fprintf(1,'Level: %d\n',lev);
    drawnow
    for coeff=1:length(cnames)
        fprintf(1,'%s...',cnames{coeff});
        drawnow
        if userefmask
            mask=wapi_readwd3(wd3mask,cnames{coeff},lev,1);
            mask=mask{1};
            refmaskidx=abs(mask)>=mean(abs(mask(mask~=0)));
            auc=zeros(1,size(mask,2));
            clear mask
        else
            auc=[];
        end
        [first,last]=wapi_idxcoef3(cnames{coeff},s,lev);
        dl=last-first+1;
        if stationary
            chunks=wapi_defaults('calculationChunks'); 
        else
            chunks=1;
        end
        if userefmask
            for p=1:chunks
                f=(p-1).*floor(dl/chunks) + 1;
                l=p .*floor(dl/chunks);
                if p==chunks
                    l=dl;
                end
                d=wapi_readwd3(wd3name,first+f-1,first+l-1,numfr);
                auc(1,f:l)=trapz(petTime,d);
                clear d
            end
            aucthr=mean(abs(auc(refmaskidx))) .* (wapi_defaults('aucThresholdPercent')/100); % aucthr=0;
            clear refmaskidx
            clear auc
        else
            aucthr=eps; % minimum computationally possible (even less actually)
        end
        fprintf(1,'AUC THRESHOLD: %.2f  ',aucthr);
        drawnow
        for p=chunks:-1:1
            f=(p-1).*floor(dl/chunks) + 1;
            l=p .*floor(dl/chunks);
            if p==chunks
                l=dl;
            end
            if chunks>1
                fprintf(1,'\nCalculating on chunk %d of %d.\n', chunks-p+1, chunks);
                drawnow
            end
            d=wapi_readwd3(wd3name,first+f-1,first+l-1,numfr);
            if interpolate
                D=interp1(petTime,d,time);
            else
                D=d;
            end
            clear d
            [tmp,tmp2]=piw_mlinlogan(D,timeL,inputL,numpoints,weights,0,k2ref,aucthr); % weights:ones, k2 reference: 0.1
            par=[tmp2;tmp];
            clear tmp tmp2
            drawnow;
            clear D
            if lev==1 && coeff==1 && p==chunks
                numP=size(par,1)-1;

                for i=1:numP
                    Wfname{i}=strcat(outfname(1:odp-1),num2str(ncoeff),'-',num2str(depth),'_mllog_p',num2str(i),'.wd3'); 
                    saveparwav(Wfname{i},lc,s);
                    wapi_extendFile([Wfname{i} '1'],lc*8);
                end
                Wfname_std=strcat(outfname(1:odp-1),num2str(ncoeff),'-',num2str(depth),'_mllog_std.wd3');
                saveparwav(Wfname_std,lc,s);
                wapi_extendFile([Wfname_std '1'],lc*8);
            end
            for i=1:numP
                fid=fopen([Wfname{i} '1'],'r+','l');
                if fid<0
                    retries=60;
                    while retries
                        pause(1);
                        fid=fopen([Wfname{i} '1'],'r+','l');
                        if fid>0
                            break;
                        end
                        retries=retries-1;
                    end
                end
                if fid<0
                    fid=fid; %#ok<ASGSL>
                end
                fseek(fid,(first+f-1)*8,-1);
                fwrite(fid,par(1+i,:),'double');

                fclose(fid);
            end
            fid=fopen([Wfname_std '1'],'r+','l');
            if fid<0
                retries=60;
                while retries
                    pause(1);
                    fid=fopen([Wfname_std '1'],'r+','l');
                    if fid>0
                        break;
                    end
                    retries=retries-1;
                end
            end
            if fid<0
                fid=fid; %#ok<ASGSL>
            end
            fseek(fid,(first+f-1)*8,-1); % only 1 is substracted since first double is to be skipped in file
            fwrite(fid,par(1,:),'double');
            fclose(fid);
        end
    end
    clear par
    fprintf(1,'done.\n');
    drawnow
end
fprintf(1,'Approximation...');
drawnow
[first,last]=wapi_idxcoef3('app',s,1);  %#ok<ASGLU>
% d=readwd3(wd3name,'app',1,numfr);
aucthr=10;%eps;
for p=chunks:-1:1
    f=(p-1).*floor(dl/chunks) + 1;
    l=p .*floor(dl/chunks);
    if p==chunks
        l=dl;
    end
    if chunks>1
        fprintf(1,'\nCalculating on chunk %d of %d.\n', chunks-p+1, chunks);
        drawnow
    end
    d=wapi_readwd3(wd3name,first+f-1,first+l-1,numfr);
    if interpolate
        D=interp1(petTime,d,time);
    else
        D=d;
    end
    clear d
    [tmp,tmp2]=piw_mlinlogan(D,timeL,inputL,numpoints,weights,0,k2ref,aucthr); % weights:ones, k2 reference: 0.1
    par=[tmp2;tmp];
    clear tmp tmp2
    clear D
    drawnow
    for i=1:numP
        fid=fopen([Wfname{i} '1'],'r+','l');
        if fid<0
            retries=60;
            while retries
                pause(1);
                fid=fopen([Wfname{i} '1'],'r+','l');
                if fid>0
                    break;
                end
                retries=retries-1;
            end
        end
        if fid<0
            fid=fid; %#ok<ASGSL> %just to be able to set a breakpoint for debugging 
        end
        fseek(fid,(first+f-1)*8,-1);
        fwrite(fid,par(1+i,:),'double');
        fclose(fid);
    end
    fid=fopen([Wfname_std '1'],'r+','l');
    if fid<0
        retries=60;
        while retries
            pause(1);
            fid=fopen([Wfname_std '1'],'r+','l');
            if fid>0
                break;
            end
            retries=retries-1;
        end
    end
    if fid<0
        fid=fid; %#ok<ASGSL>
    end
    fseek(fid,(first+f-1)*8,-1); % only 1 is substracted since first double is to be skipped in file
    fwrite(fid,par(1,:),'double');
    fclose(fid);
end
fprintf(1,'done.\n');
drawnow
clear cm

outfnames={};
for p=1:numP
    fprintf(1,'Inverting result back to normal space...');
    drawnow
    pim.volume=wapi_waverec3st(Wfname{p},[],rh,rg,rhz,rgz);
    if ~isempty(lims)
      tmp=zeros(invol.orisize);
      tmp(lims(1,1):lims(2,1),lims(1,2):lims(2,2),lims(1,3):lims(2,3))=pim.volume;
      pim.volume=tmp;
    end
    if numP>1
      pim.filename=strcat(outfname(1:odp-1),num2str(ncoeff),'-',num2str(depth),'_mllog_',sprintf('%d',p),outfname(odp:end));
    else
      pim.filename=strcat(outfname(1:odp-1),num2str(ncoeff),'-',num2str(depth),'_mllog',outfname(odp:end));
    end;

    dtype=wapi_defaults('outputDataType'); % 'int32';
    dato     = file_array(pim.filename,invol.orisize,[spm_type(dtype) spm_platform('bigend')],0,max(pim.volume(:))./spm_type(dtype,'maxval'),0);
    N        = nifti;
    N.dat    = dato;
    N.mat    = invol.mat;
    N.mat0   = invol.mat;
    N.mat_intent  = 'aligned';
    N.mat0_intent = 'scanner';
    N.descrip     = 'WAPI parametric image of DVR';

    create(N);
    fprintf(1,'Saving paramteric image in: ''%s''.\n',pim.filename);
    drawnow
    
    N.dat(:)=pim.volume(:);
    
    outfnames{end+1}=pim.filename; %#ok<AGROW>
end
if ~nargout
    clear outfnames
end
fprintf(1,'Deleting wavelet transform file of input!...');
drawnow
delete([wd3name '*']);
if userefmask
    delete([wd3mask '*']);
end
fprintf(1,'done\n');
drawnow
if deleteoutputwd3
	sstr={'','s'};
	sstr=sstr{double(numP>1)+1};
	fprintf(1,'Deleting wavelet transform file%s of output%s!...', sstr, sstr);
	drawnow
	for p=1:numP
	  delete([Wfname{p} '*']);
	end
	fprintf(1,'done\n');
	drawnow
end

return;

function saveparwav(fname,lc,s)

wapi_writevol(fname,s);
fid=fopen(fname,'a','l');
if fid<0
    retries=60;
    while retries
        pause(1);
        fid=fopen(fname,'a','l');
        if fid>0
            break;
        end
        retries=retries-1;
    end
end
if (fid<2)
  error('Error opening output file.');
end
fwrite(fid,1,'double');
fwrite(fid,lc,'double');
fclose(fid);
fid=fopen([fname '1'],'w','l');
if fid<0
    retries=60;
    while retries
        pause(1);
        fid=fopen([fname '1'],'w','l');
        if fid>0
            break;
        end
        retries=retries-1;
    end
end
fwrite(fid,1,'double');
fclose(fid);

function D=unpack_wd3(d,f,l) %#ok<DEFNU>

for fr=length(d):-1:1
    D(fr,:)=d{fr}(1,f:l); 
    if(any(isnan(D(fr,:))))
        warning('NaNs detected in coefficients (frame %d)!!',fr); %#ok<WNTAG>
    end
end

function se = strel_bol(r)

% From SPM 12
% STREL_BOL constructs a 3D sphere with the specified radius
% that can be used as structural element in 3D image processing
%
% See STREL, IMERODE, IMDILATE (image processing toolbox)

dim = [2*r+1, 2*r+1, 2*r+1];
se  = zeros(dim);
for i=1:dim(1)
  for j=1:dim(2)
    for k=1:dim(3)
      x = i-1-r;
      y = j-1-r;
      z = k-1-r;
      if norm([x y z])<=r
        se(i,j,k) = 1;
      end
    end
  end
end

