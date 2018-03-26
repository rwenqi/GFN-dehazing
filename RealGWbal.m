%º¯Êý¶¨Òå
function y=RealGWbal(Image) 
Image=double(Image);
r=Image(:,:,1); 
g=Image(:,:,2); 
b=Image(:,:,3); 
[x,y]=size(r);
avgR = mean(r(:)); 
avgG = mean(g(:)); 
avgB = mean(b(:)); 
avgRGB = [avgR avgG avgB]; 
grayValue = (avgR + avgG + avgB)/3 ;
scaleValue = grayValue./(avgRGB+0.001); 
R = scaleValue(1) * r;
G = scaleValue(2) * g;
B = scaleValue(3) * b;
for i=1:x
    for j=1:y
        if(R(i,j)>255)
            R(i,j)=255;
        end
        if(G(i,j)>255)
            G(i,j)=255;
        end
        if(B(i,j)>255)
            B(i,j)=255;
        end
    end
end
% maxB=max(B(:));
newI(:,:,1)=R;
newI(:,:,2)=G;
newI(:,:,3)=B;
y=newI; 