function block_sim_feature = block_sim_feature(image)

res_num = 5; %original image + 4 downsampled images, which means 4 self-similarities.
rb = 2^res_num;
cb = 2^res_num;
    
for tt=1:res_num
    
    [m n z] = size(image);
    if z > 1
        image = rgb2gray(image);
    end
    image=double(image);
    
    % maximum block indices
    r = floor(m/rb);
    c = floor(n/cb);
    
    signal_block=zeros(rb*cb,r*c);
    var_block=zeros(1,r*c);
    
    k=0;
    for i=1:r
        for j=1:c
            k=k+1;
            row = (rb*(i-1)+1):rb*i;
            col = (cb*(j-1)+1):cb*j;
            image_temp = image(row,col);
            
            var_block(:,k)=std(image_temp(:));
        end
    end
    
    if tt==1
        feat(tt,:) =var_block;
    else
        temp = feat(tt-1,:);
        feat(tt,:) =var_block(:,1:length(temp));
    end


    
    image = imresize(image,0.5,'bicubic');
    
    rb = rb/2;
    cb = cb/2;
end

for i=1:res_num-1
    block_sim_feature(1,i) = sum(sqrt(abs(feat(1,:)-feat(i+1,:))))/length(feat(i+1,:));

end

end

