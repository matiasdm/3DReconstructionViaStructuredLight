function [I] = SimulateProjection(P,D);

[m,n,c] = size(P);

I = P; %ini
for i = 1:m;
    for j = 1:n;
        j2 = j - D(i,j)/sqrt(2);
        i2 = i - D(i,j)/sqrt(2);
        
        % check we are inside the image, 
        if floor(i2)>0 && ceil(i2)<m+1 && floor(j2)>0 && ceil(j2)<n+1;
            k1 = (ceil(i2)-i2)*(ceil(j2)-j2) ;
            k2 = (ceil(i2)-i2)*(j2-floor(j2));
            k3 = (i2-floor(i2))*(ceil(j2)-j2);
            k4 = (i2-floor(i2))*(j2-floor(j2));
        
            if k1==0 && k2==0 && k3==0 && k4==0;
                I(i,j,:) = P(uint8(i2),uint8(j2),:);
            else
                I(i,j,:) = ( k4 * P(ceil(i2),ceil(j2),:  ) ...
                         +   k3 * P(ceil(i2),floor(j2),: ) ...
                         +   k2 * P(floor(i2),ceil(j2),: ) ...
                         +   k1 * P(floor(i2),floor(j2),:  ) ) ...
                         / ( k1+k2+k3+k4 );
            end
        else 
            I(i,j,:) = 0;
        end
    end
end

figure; imshow(P,[]);
figure; imshow(uint8(I),[]);