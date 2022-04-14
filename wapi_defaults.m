function val=wapi_defaults(name)
%wapi_defaults          Default values for adjustable parameters
%
% val=wapi_defaults(name)
%
% Get the default values for adjustable parameters used during the WAPI
% calculations. You can edit this file to change the defaults.
%
% Input:
%
% name          Name of default parameter. Please <a href="matlab: type wapi_defaults">look at the contents of</a>
%               <a href="matlab: type wapi_defaults">this .m file</a> to see the list of valid values and explanation.
%
%
% Output:
%
% val           The value for the required parameter.
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


val='';
if ~nargin || ~ischar(name)
    return;
end
switch name
    case 'outputDataType'
        % This will determine the numeric representation of data in the
        % output file. The values in the output WAPI image are scaled to
        % use the maximum range of the specified data type. Any valid SPM
        % data type string is valid here.
        val='int32';
    case 'aucThresholdPercent'
        % This determines which coefficients are excluded from Logan
        % fitting as described in the help for wapi_wapi (see part on input
        % argument refmaskfile).
        % Tha value should be numerical % value, i.e. in the range of
        % 0 to 100.
        val=1;
    case 'calculationChunks'
        % This determines in how many chunks the stationary wavelet
        % transformed PET data is processed when calculating the Logan fit
        % to the coefficient curves. In case of dyadic WT there is always
        % only a single chunk.
        % Thus you can adjust the granularity of calculations. Set to lower
        % value if you have lot of memory especially on 64-bit systems.
        val=1; %30;
end
