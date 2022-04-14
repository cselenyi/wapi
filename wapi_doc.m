function wapi_doc
%
% 3rd party uility functions used by WAPI toolbox:
% <a href="matlab: doc piw_mlinlogan">piw_mlinlogan</a>          Multi-linear (reference-)Logan plot (PIWAVE 3.1)
% <a href="matlab: doc piw_mlinlogan">piw_mlinlogan_mex</a>      Multi-linear Logan computational core (PIWAVE 3.1)
% <a href="matlab: doc uvi_lemarie">uvi_lemarie</a>            Get quadrature filters given by Battle and Lemarie (Uvi_Wave.300)
%
% Proper WAPI toolbox functions (toolbox manifest):
% <a href="matlab: doc wapi_defaults">wapi_defaults</a>          Default values for adjustable parameters.
% <a href="matlab: doc wapi_dwt3dyn">wapi_dwt3dyn</a>           Wavelet transform 4D PET image frame-by-frame
% <a href="matlab: doc wapi_extendFile">wapi_extendFile</a>        Extend file to specified size
% <a href="matlab: doc wapi_extendFile">wapi_extendFile_mex</a>    Processing core for extending file
% <a href="matlab: doc wapi_getframe">wapi_getframe</a>          Read 3D image data for a single frame
% <a href="matlab: doc wapi_gui">wapi_gui</a>               Graphical interface for WAPI
% <a href="matlab: doc wapi_waverec3st">wapi_idwt3mex_v7</a>       Isotrop 3D inverse discrete wavelet transform
% <a href="matlab: doc wapi_waverec3st">wapi_idwt3z_xhigh_mex</a>  Anisotrop 3D inverse discrete wavelet transform
%                        (only high-pass filtering along X-axis)
% <a href="matlab: doc wapi_waverec3st">wapi_idwt3z_xlow_mex</a>   Anisotrop 3D inverse discrete wavelet transform
%                        (only low-pass filtering along X-axis)
% <a href="matlab: doc wapi_waverec3st">wapi_idwt3zmex_v7</a>      Anisotrop 3D inverse discrete wavelet transform
% <a href="matlab: doc wapi_idxcoef3">wapi_idxcoef3</a>          Extract 3D approximation/detail coefficients
% <a href="matlab: doc wapi_matchIF2Frames">wapi_matchIF2Frames</a>    Match input function to PET frame timing
% <a href="matlab: doc wapi_waverec3st">wapi_mexFillArrayFromFile</a> Fill array data from file (used by <a href="matlab: doc wapi_waverec3st">wapi_waverec3st</a>)
% <a href="matlab: doc wapi_wavedec3st">wapi_mthread_dwt3mex</a>   Multi-threaded isotrop 3D forward discrete wavelet
%                        transform
% <a href="matlab: doc wapi_wavedec3st">wapi_mthread_dwt3zmex</a>  Multi-threaded anisotrop 3D forward discrete
%                        wavelet transform
% <a href="matlab: doc wapi_readblood">wapi_readblood</a>         Read blood input file.
% <a href="matlab: doc wapi_readref">wapi_readref</a>           Read reference region radioactivity file.
% <a href="matlab: doc wapi_readvol">wapi_readvol</a>           Read volume from simple ND data file
% <a href="matlab: doc wapi_readwd3">wapi_readwd3</a>           Read header or contents from .wd3 file
% <a href="matlab: doc wapi_wapi">wapi_wapi</a>              Wavelet Aided Parametric Imaging (WAPI)
% <a href="matlab: doc wapi_wavedec3st">wapi_wavedec3st</a>        Multi-level 3D wavelet decomposition
% <a href="matlab: doc wapi_waverec3st">wapi_waverec3st</a>        Multi-level 3D wavelet reconstruction
% <a href="matlab: doc wapi_writevol">wapi_writevol</a>          Write ND volume to disk in simple 3D data file
% <a href="matlab: doc wapi_wt_trsize">wapi_wt_trsize</a>         Size of wavelet filtered data
% <a href="matlab: doc wapi_readim">wapi_readim</a>            Read 3D/4D (PET) image info.
% <a href="matlab: doc wapi_readim_nifti">wapi_readim_nifti</a>      Read 3D/4D NIfTI-1 image info

% Output:
%   3rd party uility functions used by WAPI toolbox:
%   piw_mlinlogan          Multi-linear (reference-)Logan plot (PIWAVE 3.1)
%   piw_mlinlogan_mex      Multi-linear Logan computational core (PIWAVE 3.1)
%   uvi_lemarie            Get quadrature filters given by Battle and Lemarie (Uvi_Wave.300)
%  
%   Proper WAPI toolbox functions (toolbox manifest):
%   wapi_defaults          Default values for adjustable parameters
%   wapi_dwt3dyn           Wavelet transform 4D PET image frame-by-frame
%   wapi_extendFile        Extend file to specified size
%   wapi_extendFile_mex    Processing core for extending file
%   wapi_getframe          Read 3D image data for a single frame
%   wapi_gui               Graphical interface for WAPI
%   wapi_idwt3mex_v7       Isotrop 3D inverse discrete wavelet transform
%   wapi_idwt3z_xhigh_mex  Anisotrop 3D inverse discrete wavelet transform
%                          (only high-pass filtering along X-axis)
%   wapi_idwt3z_xlow_mex   Anisotrop 3D inverse discrete wavelet transform
%                          (only low-pass filtering along X-axis)
%   wapi_idwt3zmex_v7      Anisotrop 3D inverse discrete wavelet transform
%   wapi_idxcoef3          Extract 3D approximation/detail coefficients
%   wapi_matchIF2Frames    Match input function to PET frame timing
%   wapi_mexFillArrayFromFile Fill array data from file (used by wapi_waverec3st)
%   wapi_mthread_dwt3mex   Multi-threaded isotrop 3D forward discrete wavelet
%                          transform
%   wapi_mthread_dwt3zmex  Multi-threaded anisotrop 3D forward discrete
%                          wavelet transform
%   wapi_readblood         Read blood input file.
%   wapi_readref           Read reference region radioactivity file.
%   wapi_readvol           Read volume from simple ND data file
%   wapi_readwd3           Read header or contents from .wd3 file
%   wapi_wapi              Wavelet Aided Parametric Imaging (WAPI)
%   wapi_wavedec3st        Multi-level 3D wavelet decomposition
%   wapi_waverec3st        Multi-level 3D wavelet reconstruction
%   wapi_writevol          Write ND volume to disk in simple 3D data file
%   wapi_wt_trsize         Size of wavelet filtered data
%   wapi_readim            Read 3D/4D (PET) image info.
%   wapi_readim_nifti      Read 3D/4D NIfTI-1 image info

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

doc('wapi_doc');
