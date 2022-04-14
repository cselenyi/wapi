function outIF=wapi_matchIF2Frames(IF, time)
%wapi_matchIF2Frames    Match input function to PET frame timing
%
% outIF=wapi_matchIF2Frames(IF, time) matches and linearly resamples the
% input function to the frame timing of a PET image. However, it keeps the
% time corresponding to the peak activity in the input function.
%
% Inputs:
%
% IF            The input function to be resampled. It is an Nx2 matrix
%               with the first column corresponding to time and the second
%               to (radio)activity.
%
% time          The target timing. It should have the same units as
%               IF(:,1).
%
% Output:
%
% outIF         The resampled input function with outIF(:,1) being the new
%               input function timing (equals to time except that the peak
%               time is kept) and outIF(:,2) the resampled activity.
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

inT=IF(:,1);
inV=IF(:,2);

outT=time;

[mx,idx]=max(inV);

if inT(idx)<outT(1)
    warning('Time of input function peak is before first frame-time to match!\nPeak will be ignored.'); %#ok<WNTAG>
else
    outT=[outT(outT<inT(idx)); inT(idx); outT(outT>inT(idx))];
end

if outT(1)<inT(1)
    inT=[outT(1);inT];
    inV=[inV(1);inV];
end
if outT(end)>inT(end)
    inT=[inT;outT(end)];
    inV=[inV;inV(end)];
end

outV=interp1(inT,inV,outT);

outIF=[outT outV];

return