%=====================================================================
% File: marziliano.m
% Original code written by Rony Ferzli, IVU Lab (http://ivulab.asu.edu)
% Last Revised: November 2010 by Lina Karam
%===================================================================== 
% Copyright Notice:
%
% Copyright (c) 2009-2010 Arizona Board of Regents. 
% All Rights Reserved.
%
% Contact: Lina Karam (karam@asu.edu) and Adithya V. Murthy (adithya.murthy@asu.edu) 
% 
% Image, Video, and Usabilty (IVU) Lab, ivulab.asu.edu
% Arizona State University
%
% This copyright statement may not be removed from this file or from 
% modifications to this file.
% This copyright notice must also be included in any file or product 
% that is derived from this source file. 
% 
% Redistribution and use of this code in source and binary forms, 
% with or without modification, are permitted provided that the 
% following conditions are met:  
%
% - Redistribution's of source code must retain the above copyright 
% notice, this list of conditions and the following disclaimer. 
%
% - Redistribution's in binary form must reproduce the above copyright 
% notice, this list of conditions and the following disclaimer in the 
% documentation and/or other materials provided with the distribution. 
%
% - The Image, Video, and Usability Laboratory (IVU Lab, 
% http://ivulab.asu.edu) is acknowledged in any publication that 
% reports research results using this code, copies of this code, or 
% modifications of this code.  
%
% The code and our papers are to be cited in the bibliography as:
%
% A. V. Murthy and L. J. Karam, "VBQUEST- Visual Blur QUality Evaluation SofTware", 
% http://ivulab.asu.edu/Quality/VBQUEST
%
% A. V. Murthy and L. J. Karam, "MATLAB Based Framework For Image and Video Quality Evaluation,"
% International Workshop on Quality of Multimedia Experience (QoMEX), pages 242-247, June 2010.
%
% A. V. Murthy, "MATLAB Based Framework For Image and Video Quality Evaluation," Master's thesis, Arizona State University, May 2010.
%
% DISCLAIMER:
%
% This software is provided by the copyright holders and contributors 
% "as is" and any express or implied warranties, including, but not 
% limited to, the implied warranties of merchantability and fitness for 
% a particular purpose are disclaimed. In no event shall the Arizona 
% Board of Regents, Arizona State University, IVU Lab members, or 
% contributors be liable for any direct, indirect, incidental, special,
% exemplary, or consequential damages (including, but not limited to, 
% procurement of substitute goods or services; loss of use, data, or 
% profits; or business interruption) however caused and on any theory 
% of liability, whether in contract, strict liability, or tort 
% (including negligence or otherwise) arising in any way out of the use 
% of this software, even if advised of the possibility of such damage. 
% 
%===================================================================== 

function metric = marziliano(A)
%Implementation of Marziliano Blurring Metric
A = double(A);
%E = edge(A,'Sobel',[],'vertical');
E = edge(A,'canny');
[Gx, Gy] = gradient(A);
% Magnitude
graA = abs(Gx) + abs(Gy);
[M N] = size(A);
 for m=1:M
     for n=1:N
            if (Gx(m,n)~=0)
                angle_A(m,n) = atan2(Gy(m,n),Gx(m,n))*(180/pi); % in degrees
            end
            if (Gx(m,n)==0 && Gy(m,n)==0)
                angle_A(m,n) = 0;
            end
            if (Gx(m,n)==0 && Gy(m,n)==pi/2)
              angle_A(m,n) = 90;
            end
     end
 end

% quantize the angle
angle_Arnd = 45*round(angle_A./45);
width_loc = 0;
count = 0;
for m=2:M-1
    for n=2:N-1
        if (E(m,n)==1)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%% If gradient = 0 %%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (angle_Arnd(m,n) ==180 || angle_A(m,n) ==-180)
                count = count + 1;
                    for k=0:100
                    posy1 = n-1 -k;
                    posy2 = n-2 -k;
                    if ( posy2<=0)
                        break;
                    end

                    if ((A(m,posy2) - A(m,posy1))<=0)
                        break;
                    end
                end
                width_count_side1 = k + 1 ;
                for k=0:100
                    negy1 = n+1 + k;
                    negy2 = n+2 + k;

                    if (negy2>N)
                        break;
                    end

                    if ((A(m,negy2) - A(m,negy1))>=0)
                        break;
                    end

                end
                width_count_side2 = k + 1 ;
                width_loc = [width_loc width_count_side1+width_count_side2];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%% If gradient = 0 %%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (angle_Arnd(m,n) ==0)
                count = count + 1;
                    for k=0:100
                    posy1 = n+1 +k;
                    posy2 = n+2 +k;
                    if ( posy2>N)
                        break;
                    end

                    if ((A(m,posy2) - A(m,posy1))<=0)
                        break;
                    end
                end
                width_count_side1 = k + 1 ;
                for k=0:100
                    negy1 = n-1 - k;
                    negy2 = n-2 - k;

                    if (negy2<=0)
                        break;
                    end

                    if ((A(m,negy2) - A(m,negy1))>=0)
                        break;
                    end

                end
                width_count_side2 = k + 1 ;
                width_loc = [width_loc width_count_side1+width_count_side2];
            end


        end
    end
end


s = sum(width_loc)/(count);
metric = s;
