function [offs,numfr,clen,s]=wapi_readwd3(fname,map,level,lastframe)
%wapi_readwd3           Read header or contents from .wd3 file
%
% [offs,numfr,clen,s]=wapi_readwd3(fname) reads the header information from
% a .wd3 file(set).
%
% [c,s]=wapi_readwd3(fname,1)
% [c,s]=wapi_readwd3(fname,1,[],lastframe)
% reads full contents of the .wd3 file(set) into memory. For multi-frame
% transforms all frames are by default read. Specifying argument lastframe
% ensures that only up to the given number of frames of the transform are
% read into memory.
%
% [d,s]=wapi_readwd3(fname,band,level)
% [d,s]=wapi_readwd3(fname,band,level,lastframe)
% reads a single coefficient subband or all details coefficients on a given
% level into memory for all frames (default) or for certain frames only
% (if lastframe option provided).
%
% [d,s]=wapi_readwd3(fname,first,last)
% [d,s]=wapi_readwd3(fname,first,last,lastframe)
% reads an arbitrary block of coefficients from the decomposition vector
% stored in the .wd3 file into memory for all frames (default) or for
% certain frames only (if lastframe option provided).
%
% Inputs:
%
% fname         The name of the .wd3 file to be read from.
%
% band          The name of the subband to be read from file: 'app', 'h_l',
%               'v_l', 'd_l', 'a_h', 'h_h', 'v_h', 'd_h'. See help for
%               <a href="matlab: help wapi_wavedec3st">wapi_wavedec3st</a> for more info. Additionaly, you can specify
%               'c' or 'compact' to read all details coefficients at a
%               given level.
%
% level         The level for which the coefficients are to be read.
%               Ignored in case of approximation coefficients (i.e. when
%               band='app').
%
% first         The first coefficient to be read for an arbitrary block.
%
% last          The last coefficient to be read for an arbitrary block.
%
% lastframe     If provided then coefficients for only this number of
%               frames is read into memory.
%
% Outputs:
%
% offs          The offset of the first coefficient in the .wd3 file(set).
%
% numfr         The number of frames stored in the .wd3 file(set).
%
% clen          The length of the decomposition structure vector c (see
%               <a href="matlab: help wapi_wavedec3st">wapi_wavedec3st</a> for more info).
%
% c             A cell array containing one vector per frame corresponding
%               to the decomposition structure vector c (see <a href="matlab: help wapi_wavedec3st">wapi_wavedec3st</a> 
%               for more info).
%
% d             A cell array containing one vector per frame corresponding
%               to requested part of the decomposition structure vector c
%               (specific subband or arbitrary block).
%
% s             The decomposition bookkeeping matrix s (see <a href="matlab: help wapi_wavedec3st">wapi_wavedec3st</a>
%               for more info.
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

if nargin<2
  map=[];
end
if nargin<3
    level=1;
end
if nargin<4
    lastframe=[];
end
[s,offs]=wapi_readvol(fname);

fid=fopen(fname,'r','l');
if fid<0
  error('Error opening file');
end
fseek(fid,offs,-1);
numfr=fread(fid,1,'double');
if nargin>3 && ~isempty(lastframe)
    if lastframe(1)<numfr
        numfr=round(lastframe(1));
    end
    if numfr<=0
        error('Last frame specified erroneously!');
    end
end    
clen=fread(fid,1,'double');
offs=offs+(2*8);
fclose(fid);
if ~isempty(map)
    if ischar(map)
        [first,last]=wapi_idxcoef3(map,s,level);
        for fr=1:numfr
            readwd3_d{fr}=zeros(1,last-first+1); %#ok<AGROW>
        end

        for fr=1:numfr
            fid2=fopen(strcat(fname,num2str(fr)),'r','l');
            filefr=fread(fid2,1,'double');
            if filefr~=fr
                warning('Incorrect frame number in file: ''%s''.', strcat(fname,num2str(fr))); %#ok<WNTAG>
            end
            fseek(fid2,(first-1).*8,0);
            readwd3_d{fr}(1,:) = fread(fid2,last-first+1,'double')'; 
            fclose(fid2);
        end
        offs=readwd3_d;
        numfr=s;
    else
        if nargin>2 % read only certain subset of coeffs
            first=map;
            last=level;
            if numfr>1
                readwd3_d=zeros(numfr,last-first+1); 
            end
            for fr=1:numfr
                fid2=fopen(strcat(fname,num2str(fr)),'r','l');
                filefr=fread(fid2,1,'double');
                if filefr~=fr
                    warning('Incorrect frame number in file: ''%s''.', strcat(fname,num2str(fr))); %#ok<WNTAG>
                end
                fseek(fid2,(first-1).*8,0);
                readwd3_d(fr,:) = fread(fid2,last-first+1,'double')'; 
                fclose(fid2);
            end
            offs=readwd3_d;
            numfr=s;
        else
            c(numfr,clen)=0; % preoccupy memory
            for fr=1:numfr
                fid2=fopen(strcat(fname,num2str(fr)),'r','l');
                filefr=fread(fid2,1,'double');
                if filefr~=fr
                    warning('Incorrect frame number in file: ''%s''.', strcat(fname,num2str(fr))); %#ok<WNTAG>
                end
                c(fr,:) = fread(fid2,clen,'double')';
                fclose(fid2);
            end
            offs=c;     % First output is 'c' in this case.
            numfr=s;    % Second output is 's' in this case.
            clen=[];    % Third output is undefined.
            s=[];       % Fourth output is undefined.
        end
    end
end
return;
