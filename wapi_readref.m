function ref=wapi_readref(refin)
%wapi_readref           Read reference region radioactivity file.
%
% ref=wapi_readref(filename)
%
% Read text file containing radioactivity in (reference) region.
%
% Input file should contain radioactivity only. No headerlines are assumed
% and a single column of numbers. 
%
% Input:
%
% filename      Name of file to be read.
%
%
% Output:
%
% ref           Vector of radioactivity
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


ref=textread(refin,'%f','headerlines',0,'delimiter','\t'); %#ok<DTXTRD>

return;
