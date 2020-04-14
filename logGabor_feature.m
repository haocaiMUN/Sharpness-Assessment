function gf  = logGabor_feature(img_input)

  gf = logGabor(img_input);
  
function   D= logGabor(im)
num_scale =3 ;   
num_orien =10;   
minWaveLength =3;   
mult = 2;       
sigmaOnf = 0.65;     
dThetaOnSigma =1.5;  

gf = gaborconvolve(im,  num_scale, num_orien, minWaveLength, mult, sigmaOnf, dThetaOnSigma);
%%
A=[];
C=[];

for i = 1:num_scale    
    for j = 1:num_orien    
           
         sub_im = log(abs(gf{i,j})+.1); 
           kk=abs(gf{i,j}).^2;
           C  = [C mean(kk(:))];
           
           win = [1 0 -1; 0 0 0; -1 0 1;];
           deri_im_7 = filter2(win, sub_im, 'same');
           deri_im_7=deri_im_7(:);
        
        [sigma_7 alpha_7] = gaussian_para_esti(deri_im_7);
        A=[A sigma_7 alpha_7];
    end
end
 D=[A  C];
