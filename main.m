I=imread('image.png');
M=imread('mask.png');
M1=findVesselLumen(I,M,12);
M2=M;
M2(M1>0)=0;
figure, imshow([I,M,M2],[])
B=labeloverlay(I,M2);
figure, imshow(B)