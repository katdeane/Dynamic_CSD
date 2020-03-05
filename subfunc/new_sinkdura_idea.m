load('Anesthetizedpre_Data.mat')
rawCSD= Data.GXL03.CSD{4};
CSDleaf= Data.GXL03.CSD{4};
figure, imagesc(CSDleaf),colormap(jet)
th= -0.0001; %discuss threshold level with max vs. -0.0002 (captures more signal)

CSDleaf(CSDleaf >= th)=0;
CSDleaf(:,1:200)=0;
CSDleaf(CSDleaf <= th)=1;

CSDleaf = im2bw(CSDleaf);
CSDleaf = bwareaopen(CSDleaf,50);
figure, imagesc(CSDleaf), colormap(jet)


stats=regionprops(CSDleaf);

