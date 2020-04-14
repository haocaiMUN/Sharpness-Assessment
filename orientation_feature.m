function orientation_hist = orientation_feature( img0, nnum )

np = nnum; % neighbor pixels
otr = 6; % threshold for the similar preferred orientation

kx = [0 0; -1 1];
ky = kx';

sps=zeros(np,2);
as = 2*pi/np;
r = 1;
for i = 1:np
    sps(i,1) = -r*sin((i-1)*as);
    sps(i,2) = r*cos((i-1)*as);
end
bin_num = np+1;
img = imresize( img0, 1/2 );
imgd = padarray( img, [r r], 'symmetric' );
[ row, col ] = size( imgd );
    

[SZx SZy] = size(imgd);

Tx = [0 0; -1 1];
Gx = conv2(imgd, Tx, 'same');
Gy = conv2(imgd, Tx', 'same');
Gx(:,1) = imgd(:,2) - imgd(:,1);
Gx(:,SZy) = imgd(:,SZy) - imgd(:,(SZy-1));
Gy(1,:) = imgd(2,:) - imgd(1,:);
Gy(SZx,:) = imgd(SZx,:) - imgd((SZx-1),:);
Cimg=(abs(Gx) + abs(Gy))/2;
    

Cvimg = zeros( row,col );
thre = 10;
Cvimg( Cimg>=thre ) = 1;
Oimg = atan2(Gy,Gx)/pi*180;
Cvimgc = Cvimg( r+1:row-r,r+1:col-r );
Oimgc = Oimg( r+1:row-r,r+1:col-r );
    
ssr_val = zeros( row-2*r, col-2*r );
for i = 1:np
  dx = round( r+sps(i,1) );
  dy = round( r+sps(i,2) );
  Cvimgn = Cvimg( dx+1:row-2*r+dx, dy+1:col-2*r+dy );
  Oimgn = Oimg( dx+1:row-2*r+dx, dy+1:col-2*r+dy );
  Odif = abs( Oimgn - Oimgc );
  Odif( Odif>180 ) = 360 - Odif( Odif>180 );
  Odif( Cvimgc==0 ) = 360;
  last_state = zeros( row-2*r, col-2*r );
  last_state( Cvimgn+Cvimgc==0 ) = 1;
  last_state( Odif<=otr ) = 1;
  ssr_val = ssr_val + last_state;
end

% mapping
var_map = func_std_var( img, r );
orientation_hist = zeros( bin_num, 1 );
for i = 1 : bin_num
    orientation_hist(i) = sum( sum( var_map( ssr_val==i-1) ) );
end

return;
