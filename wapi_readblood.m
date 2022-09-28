function blo=wapi_readblood(filename,timeunits)
%wapi_readblood         Read blood input file.
%
% blo=wapi_readblood(filename,[timeunits])
%
% Read text file with blood (plasma) input data.
% Input file should contain time, total blood, met. corr. plasma with time
% in seconds. One headerline assumed and 3 columns of numbers. 
%
% Inputs:
%
% filename      Name of file containing blood input data.
%
% timeunits     Optional. Required time units in the output matrix. 
%               Valid choices: 'msec', 'sec', 'min' (default is 'min').
%
% Outputs:
%
% blo           Matrix of blood data with 3 columns for: time (min), total
%               blood, met. corr. plasma.
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
    timeunits='min';
end
switch timeunits
    case 'msec'
        factor=1000;
    case 'sec'
        factor=1;
    otherwise
        factor=1/60;
end

[bloodtime,blood,plasma]=textread(filename,'%f%f%f','headerlines',1,'delimiter','\t'); %#ok<DTXTRD>

blood(blood<0)=0;
plasma(plasma<0)=0;
blo=[bloodtime*factor blood plasma];

return;
