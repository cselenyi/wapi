function [first,last] = wapi_idxcoef3(o,s,n)
%wapi_idxcoef3          Extract 3D approximation/detail coefficients
%
% [first, last]=wapi_idxcoef3(o,s,n) gets the index of the detail
% coefficients (a single sub-band or all details) at a given level or of
% the approximation coefficients at the deepest level within the wavelet
% decomposition structure vector (see output c of <a href="matlab: help wapi_wavedec3st">wapi_wavedec3st</a>). 
%
% Inputs:
%
% o             The name of the subband to get the index for: 'app', 'h_l',
%               'v_l', 'd_l', 'a_h', 'h_h', 'v_h', 'd_h'. See help for
%               <a href="matlab: help wapi_wavedec3st">wapi_wavedec3st</a> for more info. Additionaly, you can specify
%               'c' or 'compact' to get the index for all details
%               coefficients at a given level.
%
% s             The decomposition bookkeeping matrix (see help for
%               <a href="matlab: help wapi_wavedec3st">wapi_wavedec3st</a> for an explanation).
%
% n             The level from which the coefficients are to be extracted.
%               The value must be an integer such that
%               1 <= n <= size(s,1)-2.
%               If the approximation part ('app') is to be obtained that
%               the value of n is ignored.
%
% Outputs:
%
% first         The index of the first requested coeffcient within the
%               decomposition vector.
%
% last          The index of the last requested coeffcient within the
%               decomposition vector.
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

if ~any((3:4)==nargin)
    error('Wrong number of inputs');
end
if ~any((1:2)==nargout)
    error('Wrong number of outputs');
end

nmax = size(s,1)-2;
if (n<1) || (n>nmax) || (n~=fix(n))
        error('Invalid level value');
end

st=0;
if all(sum(diff(s))==0)
    st=1; % trimmed transform with trimmed s matrix (cannot be zero padded!)
end
k     = size(s,1)-n;
if st>0 % trimmed st transform
    add = prod( s(end,1:3) );
    first = ( 1 + (k-2).*7 ) .* add + 1;
else
    first = s(1,1)*s(1,2)*s(1,3)+7*sum(s(2:k-1,1).*s(2:k-1,2).*s(2:k-1,3))+1;
    add   = s(k,1)*s(k,2)*s(k,3);
end
o     = lower(o);

last=0;
switch o
    case {'a_l','app'} ,
      first=1;
      if st>0
          last=add;
      else
          last=prod(s(1,:));
      end
      return;
    case 'h_l' ,
    case 'v_l' , first = first+add;
    case 'd_l' , first = first+2*add;
    case 'a_h' , first = first+3*add;
    case 'h_h' , first = first+4*add;
    case 'v_h' , first = first+5*add;
    case 'd_h' , first = first+6*add;
    case {'c','compact'}
      last = first+7*add-1;
    otherwise   
      error('Invalid argument value');
end

if ~last
    last = first+add-1;
end
return;
