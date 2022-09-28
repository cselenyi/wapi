function [vol,cnt]=wapi_readvol(fname,endi)
%wapi_readvol           Read volume from simple ND data file
%
% [vol,cnt]=wapi_readvol(fname,endi) reads a ND volume from a file
% containing data in simple volume storage format. See <a href="matlab help wapi_writevol">wapi_writevol</a> for
% details on this format.
%
% Inputs:
%
% fname         Name of the file to read from.
%
% endi          The machine data format (endiannes) of the file. See help
%               for <a href="matlab help fopen">fopen</a> for details.
%
% Outputs:
%
% vol           The ND array of data read from the file.
%
% cnt           The number of bytes read from the file.
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
    endi='l';
end
fid=fopen(fname,'r',endi);
if(fid<1)
  error('Error opening file');
end
ndim=fread(fid,1,'double');
dim=fread(fid,ndim,'double')';
vol=zeros(dim);
vol(:)=fread(fid,prod(dim),'double');
cnt=ftell(fid);
fclose(fid);
return;

