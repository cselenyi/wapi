function  [h,g,rh,rg]=uvi_lemarie(num_coefs)
%uvi_lemarie           Get quadrature filters given by Battle and Lemarie
%
% [h,g,rh,rg] = uvi_lemarie (num_coefs) returns the coefficients of
% orthonormal  Battle-Lemarie wavelets. 
%
% Inputs:
%
% num_coefs     specifies the desired number of coefficients (filter kernel
%               length).
%
% Outputs:
%
% h             decomposition low-pass filter
%
% g             decomposition high-pass filter
%
% rh            reconstruction low-pass filter
%
% rg            reconstruction high-pass filter
%
% References:
%
% S. Mallat, "A Theory for Multiresolution Signal Decomposition: The
% Wavelet Representation", IEEE Trans. on Patt. An. and Mach. Intell., July
% 1989 

%
% Based on LEMARIE and RH2RG from Uvi_Wave package.
% Modified by Zsolt Cselényi. 2022-04-14
%
%--------------------------------------------------------
% Copyright (C) 1994, 1995, 1996, by Universidad de Vigo 
%                                                      
%                                                      
% Uvi_Wave is free software; you can redistribute it and/or modify it      
% under the terms of the GNU General Public License as published by the    
% Free Software Foundation; either version 2, or (at your option) any      
% later version.                                                           
%                                                                          
% Uvi_Wave is distributed in the hope that it will be useful, but WITHOUT  
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or    
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License    
% for more details.                                                        
%                                                                          
% You should have received a copy of the GNU General Public License        
% along with Uvi_Wave; see the file COPYING.  If not, write to the Free    
% Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.             
%                                                                          
%       Author: Jose Martin Garcia
%       e-mail: Uvi_Wave@tsc.uvigo.es
%--------------------------------------------------------


% frequency axis
L=128;
w=[0:2*pi/L:2*pi*(1-1/L)];
w(1)=eps;
w(65)=w(65)+1e-15;

% calculation of frequency response of analysis lowpass filter 
num=0;den=0;
K=36;
for k=-K:K,
	num=1./((w+2*pi*k).^8)+num;
	den=1./((2*w+2*pi*k).^8)+den;
end
h=sqrt(num./(2.^8*den));
h(1)=1;

% obtain the time response of lowpass filter
h=real(ifft(h,L));
h=[ h(128-floor(num_coefs/2)+1:128) h(1:ceil(num_coefs/2))];
h=sqrt(2)/sum(h)*h;

[rh,rg,h,g]=rh2rg(fliplr(h));

function [rh,rg,h,g]=rh2rg(rh)
%rh2rg    Calculates all the filters from the synthesis lowpass
%	  in the orthogonal case.
%
%	  [RH,RG,H,G]=RH2RG(RH) begins with the synthesis lowpass
%	  filter (RH) and returns the synthesis highpass filter (RG),
%	  the analysis lowpass filter (H) and the analysis highpass
%	  filter (G).
%	
%	  It is an auxiliary function for orthogonal filters design.

%--------------------------------------------------------
% Copyright (C) 1994, 1995, 1996, by Universidad de Vigo 
%                                                      
%                                                      
% Uvi_Wave is free software; you can redistribute it and/or modify it      
% under the terms of the GNU General Public License as published by the    
% Free Software Foundation; either version 2, or (at your option) any      
% later version.                                                           
%                                                                          
% Uvi_Wave is distributed in the hope that it will be useful, but WITHOUT  
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or    
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License    
% for more details.                                                        
%                                                                          
% You should have received a copy of the GNU General Public License        
% along with Uvi_Wave; see the file COPYING.  If not, write to the Free    
% Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.             
%                                                                          
%       Author: Nuria Gonzalez Prelcic
%       e-mail: Uvi_Wave@tsc.uvigo.es
%--------------------------------------------------------

% Calculate rg from rh.

for i=1:length(rh)        
	rg(i) = -(-1)^i*rh(length(rh)-i+1);
end  

% Calculate h and g

h=rh(length(rh):-1:1);
g=rg(length(rg):-1:1);
