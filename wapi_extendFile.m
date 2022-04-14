function wapi_extendFile(filename,extrasize)
%wapi_extendFile        Extend file to specified size
%
% wapi_extendFile(filename,extrasize) attempts to quickly extend the size
% of the specified file with the given number of bytes. This is done using
% systems calls to extend the valid data length of the file. Thus the file
% is truly increased without having to write the given number of bytes to
% disk. Requires privileges since the contents of the extended file will be
% undefined and may contain data from previously deleted files on the disk.
% This step may be useful for speeding up writing to files from a
% process with large-volume output such as frame-by-frame transform of 4D
% PET images.
%  Uses wapi_extendFile_mex to perform the extension. It fails silently.
%
% Inputs:
%
% filename      Name of the file to be expanded. The file cannot be open.
%
% extrasize     The required additional size of the file in bytes.
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

try
    wapi_extendFile_mex(filename,extrasize);
catch
    estr=lasterr; %#ok<LERR>
    pos=strfind(estr,'EXTENDFILE_MEX: ');
    if isempty(pos)
        'Unknown error.';
    else
        estr=estr(pos+length('EXTENDFILE_MEX: '):end);
    end
    if ~exist(filename,'file')
        d.bytes=NaN;
    else
        d=dir(filename);
    end
    warning('wapi_extendFile:errorExtending','Could not extend file %s (current size: %d) with %d bytes:\n%s',filename,d.bytes,extrasize,estr);
    % wapi_extendFile_mex(filename,0);
end
