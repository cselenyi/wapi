function wapi_writevol(fname,vol)
%wapi_writevol          Write ND volume to disk in simple 3D data file
%
% wapi_writevol(fname,vol) write an ND volume to a file in simple volume
% storage format. The specifications of the format:
% - data are stored as doubles in little endian machine format.
% - 1st number stored is the length of the size of the volume, i.e. the
% dimensionality of the volume (=N).
% - Numbers 2 to N+1: the size of the array along dimensions 1 to N.
% - Rest of the file: the contents of the array stored linearly.
%
% Inputs:
%
% fname         Name of the file to write to.
%
% vol           The 3D array of data to be writen to the file.
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

fid=fopen(fname,'w','l');
if(fid<1)
  error('Error opening file.');
end
siz=size(vol);
fwrite(fid,length(siz),'double');
fwrite(fid,siz,'double');
fwrite(fid,vol(:),'double');
fclose(fid);
return;

