function trsize=wapi_wt_trsize(dimsize, filterlen, st)
%wapi_wt_trsize         Size of wavelet filtered data
%
% trsize=wapi_wt_trsize(dimsize, filterlen, st) returns the size of the
% wavelet filtered data for a single subband in a single iteration of the
% wavelet decomposition (i.e. a single convolution [and decimation] step
% applied along each dimension). You must specify whether dyadic or
% stationary transform is to be applied.
%
% Inputs:
%
% dimsize       Size of the input data to be transformed. Typically a
%               vector with 3 elements indicating size along X, Y and Z.
%
% filterlen     The length of the wavelet filter kernel to be used.
%
% st            Value 0 indicates dyadic transform, a value of 1 indicates
%               stationary transform.
%
% Output:
%
% trsize        The size of a given subband of filtered data. trsize is a
%               vector of length equal to dimsize.
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

if st
    trsize=dimsize;
else
    div=2;
    if numel(filterlen)==1
        filterlen=repmat(filterlen,size(dimsize));
    end
    trsize=floor( (dimsize+filterlen-1)/div);
end
