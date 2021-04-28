%Gross and detail morphology.
clear
nWMC=load('cbWCM3.mat','cb');
nWME=load('cbWEM3.mat','cb');
nWFC=load('cbWCF2.mat','cb');
nWFE=load('cbWEF.mat','cb');
nBMC=load('cbBCM.mat','cb');
nBME=load('cbBEM3.mat','cb');
nBFC=load('cbBCF.mat','cb');
nBFE=load('cbBEF2.mat','cb');
%% cb WMC B6MaleControl
rIm=nWMC.cb.WT.R.M.C{1}(5);%receptor image
dIm=nWMC.cb.WT.D.M.C{1}(5);%dapi image
im=enhanced1Im(rIm,dIm,0.2,0.2,[5 0.5],[10 5],2);
axis equal
axis off
imwrite(im,'cbWMC.tif','tiff','resolution',300)
im=enhanced1ImAltAllIm(rIm,dIm,0.5,0.3,[5 0.8],[10 5],2);axis equal;axis off
exportgraphics(gca,'cbWMC_enhanced.tif','resolution',300);

%% cb BMC   BTBR Male Control
rIm=nBMC.cb.BT.R.M.C{1}(2);
dIm=nBMC.cb.BT.D.M.C{1}(2);
im=enhanced1Im(rIm,dIm,0.1,0.1,[5 0.8],[10 5],2);

axis equal
axis off
imwrite(im,'cbBMC.tif','tiff','resolution',300)
%working on this
im=enhanced1ImAltAllIm(rIm,dIm,0.4,0.4,[10 1],[10 5],2);axis equal;axis off
exportgraphics(gca,'cbBMC_enhanced.tif','resolution',300)

%% cb WFC
rIm=nWFC.cb.WT.R.F.C{4}(3);
dIm=nWFC.cb.WT.D.F.C{4}(3);
im=enhanced1Im(rIm,dIm,0.2,0.1,[5 0.8],[10 5],2);


%detail showing basket around a PC soma
%im=enhanced1ImAlt(rIm,dIm,0.2,0.2,[10 0.2],[5 2],2);axis equal

axis equal
axis off
imwrite(im,'cbWFC.tif','tiff','resolution',300)

im=enhanced1ImAltAllIm(rIm,dIm,0.45,0.4,[10 1],[10 5],2);axis equal; axis off
exportgraphics(gca,'cbWFC_enhanced.tif','resolution',300)

%detail showing basket around a PC soma
im=enhanced1ImAlt(rIm,dIm,0.2,0.2,[10 0.8],[20 5],[488 733],2);axis equal;axis off
exportgraphics(gca,'cbWFC_enhanced_detailed.tif','resolution',300)




%% cb BFC
rIm=nBFC.cb.BT.R.F.C{2}(3);
dIm=nBFC.cb.BT.D.F.C{2}(3);
im=enhanced1Im(rIm,dIm,0.1,0.1,[5 0.8],[10 5],2);

axis equal
axis off
imwrite(im,'cbBFC.tif','tiff','resolution',300)

im=enhanced1ImAltAllIm(rIm,dIm,0.1,0.4,[10 1],[10 5],2);axis equal; axis off
exportgraphics(gca,'cbBFC_enhanced.tif','resolution',300);

im=enhanced1ImAlt(rIm,dIm,0.2,0.7,[10 0.2],[20 2],[567 413],2);axis equal;axis off
exportgraphics(gca,'cbBFC_enhanced_detailed.tif','resolution',300)


%% NNMDA
clear
nWMC=load('nmdaWCM.mat','nmda');
nWME=load('nmdaWEM2.mat','nmda');%ols nmdsWEM.mat
nWFC=load('nmdaWCF2.mat','nmda');
nWFE=load('nmdaWEF2.mat','nmda');
nBMC=load('nmdaBCM.mat','nmda');
nBME=load('nmdaBEM3.mat','nmda');%old nmdaBEM2.mat
nBFC=load('nmdaBCF.mat','nmda');
nBFE=load('nmdaBEF.mat','nmda');

%% nm WMC
rIm=nWMC.nmda.WT.R.M.C{3}(3);
dIm=nWMC.nmda.WT.D.M.C{3}(3);
im=enhanced1Im(rIm,dIm,0.2,0.2,[10 0.2],[10 5],1);

axis equal
axis off
imwrite(im,'nmWMC.tif','tiff','resolution',300)

im=enhanced1ImAltAllIm(rIm,dIm,0.4,0.1,[20 1],[10 5],1);axis equal;axis off
exportgraphics(gca,'nmWMC_enhanced.tif','resolution',300);


%% nm BMC
rIm=nBMC.nmda.BT.R.M.C{3}(3);
dIm=nBMC.nmda.BT.D.M.C{3}(3);
im=enhanced1Im(rIm,dIm,0.2,0.6,[3 0.2],[10 5],1);
axis equal
axis off
imwrite(im,'nmBMC.tif','tiff','resolution',300)

im=enhanced1ImAltAllIm(rIm,dIm,0.3,0.1,[20 1],[10 5],1);axis equal;axis off
exportgraphics(gca,'nmBMC_enchanced.tif','resolution',300);

%% nm WFC
rIm=nWFC.nmda.WT.R.F.C{3}(1);
dIm=nWFC.nmda.WT.D.F.C{3}(1);
im=enhanced1Im(rIm,dIm,0.2,0.75,[10 2],[10 5],1);
axis equal
axis off
imwrite(im,'nmWFC.tif','tiff','resolution',300)

im=enhanced1ImAltAllIm(rIm,dIm,0.4,0.1,[20 1],[10 5],1);axis equal;axis off
exportgraphics(gca,'nmWFC_enhanced.tif','resolution',300);

%% nm BFC
rIm=nBFC.nmda.BT.R.F.C{4}(2);
dIm=nBFC.nmda.BT.D.F.C{4}(2);
im=enhanced1Im(rIm,dIm,0.2,0.2,[10 0.2],[10 5],1);

axis equal
axis off
imwrite(im,'nmBFC.tif','tiff','resolution',300)
im=enhanced1ImAltAllIm(rIm,dIm,0.3,0.1,[10 1],[10 5],1);axis equal;axis off
exportgraphics(gca,'nmBFC_enhanced.tif','resolution',300);

%% The z-stacks

%% cb BMC

clear
homedir='/home/fidel/Daniela/BTBR/';
databasedir='/home/Daniela.Monje/paperimages/v2_images/';
 load DetailedMorphIm
cd(homedir)
%good. Use this for Fig 2.

%% cb BMC
[R,D]=getmeIm(dirnamesBT{1}(3));
imR.Data.imbc=R(:,:,1:10);%mpR;%apR./max(apR(:));
imD.Data.imbc=(D(:,:,1:8));%max(M(:,:,3,:,a),[],4);
% [rim,dim,x,y]=getTotFIm(imR,imD,[],200,[0.1 0.999],[0.2 0.8],100,20,'g',[10 1]); 
x=689;y=325;
 [bcMCB,dim,x,y,rgbAllmax,rgbAllsum]=getTotFIm(imR,imD,[x,y],300,[0.4 0.998],[0.3 0.99],-1,30,'g',[10 1],[20 5]);
imwrite(dim,'cbBMC_3d_detail.tif','tiff','resolution',300)
imwrite(rgbAllsum,'cbBMC_3d_all.tif','tiff','resolution',300)


%% not good

 [R,D]=getmeIm(dirnamesBT{7}(4));
imR.Data.imbc=R(:,:,1:10);%mpR;%apR./max(apR(:));
imD.Data.imbc=(D(:,:,1:8));%max(M(:,:,3,:,a),[],4);
 %[rim,dim,x,y]=getTotFIm(imR,imD,[],200,[0.1 0.999],[0.2 0.8],100,20,'g',[10 1]); 
x=595;y=609;

 [bcMCB,dim,x,y]=getTotFIm(imR,imD,[x,y],300,[0.3 1],[0.5 0.7],10,60,'g',[10 3],[20 5]);

 %% cb BFC
 [R,D]=getmeIm(dirnamesBT{7}(3));
imR.Data.imbc=R(:,:,1:10);%mpR;%apR./max(apR(:));
imD.Data.imbc=(D(:,:,1:8));%max(M(:,:,3,:,a),[],4);
% [rim,dim,x,y]=getTotFIm(imR,imD,[],300,[0.1 0.999],[0.2 0.8],0,60,'g',[10 2],[20 10]); 
%x=588;y=357;
x=540;y=600;
 [bcMCB,dim,x,y,rgbAllmax,rgbAllsum]=getTotFIm(imR,imD,[x,y],300,[0.4 0.98],[0.4 0.998],0,100,'g',[10 2],[20 10]);

imwrite(dim,'cbBFC_3d_detail.tif','tiff','resolution',300)
imwrite(rgbAllsum,'cbBFC_3d_all.tif','tiff','resolution',300)



 %% nm BFC
   [R,D]=getmeIm(dirnamesBT{3}(1));
imR.Data.imbc=R(:,:,1:10);%mpR;%apR./max(apR(:));
imD.Data.imbc=(D(:,:,1:8));%max(M(:,:,3,:,a),[],4);
% [rim,dim,x,y]=getTotFIm(imR,imD,[],200,[0.1 0.999],[0.2 0.8],100,20,'r',[10 1],[20 10]); 
x=530;y=570;
 [bcMCB,dim,x,y,rgbAllmax,rgbAllsum]=getTotFIm(imR,imD,[x,y],300,[0.3 0.9],[0.3 0.99],0,100,'r',[10 3],[20 7]);
 imwrite(dim,'nmBFC_3d_detail.tif','tiff','resolution',300)
 imwrite(rgbAllsum,'nmBFC_3d_all.tif','tiff','resolution',300)
 

%%
   [R,D]=getmeIm(dirnamesBT{3}(2));
imR.Data.imbc=R(:,:,3:5);%mpR;%apR./max(apR(:));
imD.Data.imbc=(D(:,:,1:8));%max(M(:,:,3,:,a),[],4);
 [rim,dim,x,y]=getTotFIm(imR,imD,[],300,[0.1 0.999],[0.2 0.8],100,20,'r',[10 1],[20 10]); 
%x=419;y=323;
 [bcMCB,dim,x,y,rgbAllmax,rgbAllsum]=getTotFIm(imR,imD,[x,y],300,[0.2 0.99],[0.3 0.99],10,60,'r',[10 3],[20 10]);

 %% nm BMC %nice dendrite
    [R,D,xlmd]=getmeIm(dirnamesBT{6}(2));
imR.Data.imbc=R(:,:,1:10);%mpR;%apR./max(apR(:));
imD.Data.imbc=(D(:,:,1:8));%max(M(:,:,3,:,a),[],4);
 %[rim,dim,x,y]=getTotFIm(imR,imD,[],200,[0.1 0.999],[0.2 0.8],100,20,'r',[10 1]); 
x=583;y=540;
 [bcMCB,dim,x,y,rgbAllmax,rgbAllsum]=getTotFIm(imR,imD,[x,y],300,[0.4 0.999],[0.2 0.99],0,50,'r',[10 2],[20 7]);
 imwrite(dim,'nmBMC_3d_detail.tif','tiff','resolution',300)
  imwrite(rgbAllsum,'nmBMC_3d_all.tif','tiff','resolution',300)


 %% nm WMC %showming main dendrites
[R,D]=getmeIm(dirnamesWT{4}(2));
imR.Data.imbc=R(:,:,1:10);%mpR;%apR./max(apR(:));
imD.Data.imbc=(D(:,:,1:8));%max(M(:,:,3,:,a),[],4);
% [rim,dim,x,y]=getTotFIm(imR,imD,[],200,[0.1 0.999],[0.2 0.8],100,20,'r',[10 2]); 
x=476;y=365;
 [bcMCB,dim,x,y,rgbAllmax,rgbAllsum]=getTotFIm(imR,imD,[x,y],300,[0.4 0.99],[0.2 0.99],0,100,'r',[10 2],[20 7]);
 imwrite(dim,'nmWMC_3d_detail.tif','tiff','resolution',300)
 imwrite(rgbAllsum,'nmWMC_3d_all.tif','tiff','resolution',300)


 %% nm WFC
 [R,D]=getmeIm(dirnamesWT{3}(1)); 
imR.Data.imbc=R(:,:,1:10);%mpR;%apR./max(apR(:));
imD.Data.imbc=(D(:,:,1:5));%max(M(:,:,3,:,a),[],4);
 %[rim,dim,x,y]=getTotFIm(imR,imD,[],200,[0.1 0.999],[0.2 0.8],100,20,'r',[10 1],[20 10]); 
x=518;y=587;
 [bcMCB,dim,x,y,rgbAllmax,rgbAllsum]=getTotFIm(imR,imD,[x,y],300,[0.4 0.99],[0.3 0.99],0,100,'r',[10 2],[20 5]);
imwrite(dim,'nmWFC_3d_detail.tif','tiff','resolution',300)
imwrite(rgbAllsum,'nmWFC_3d_all.tif','tiff','resolution',300)

%% cb WMC %super good
  [R,D]=getmeIm(dirnamesWT{7}(3));
imR.Data.imbc=R(:,:,4:8);%mpR;%apR./max(apR(:));
imD.Data.imbc=(D(:,:,1:10));%max(M(:,:,3,:,a),[],4);
 %[rim,dim,x,y]=getTotFIm(imR,imD,[],300,[0.1 0.999],[0.2 0.8],100,20,'g',[10 1],[20 10]); 
x=663;y=327;
 [bcMCB,dim,x,y,rgbAllmax,rgbAllsum]=getTotFIm(imR,imD,[x,y],300,[0.5 0.99],[0.4 0.999],0,50,'g',[20 1],[20 5]);
 imwrite(dim,'cbWMC_3d_detail.tif','tiff','resolution',300)
 imwrite(rgbAllsum,'cbWMC_3d_all.tif','tiff','resolution',300)

 %% cb WFC
  [R,D]=getmeIm(dirnamesWT{5}(2));
imR.Data.imbc=R(:,:,1:10);%mpR;%apR./max(apR(:));
imD.Data.imbc=(D(:,:,1:10));%max(M(:,:,3,:,a),[],4);
 %[rim,dim,x,y]=getTotFIm(imR,imD,[],200,[0.6 0.999],[0.2 0.8],0,30,'g',[10 2],[20 10]); 
x=439;y=724;
x=527;y=680
 [bcMCB,dim,x,y,rgbAllmax,rgbAllsum]=getTotFIm(imR,imD,[x,y],300,[0.3 0.998],[0.05 0.99],0,50,'g',[10 2],[20 5]);
imwrite(dim,'cbWFC_3d_detail.tif','tiff','resolution',300)
imwrite(rgbAllsum,'cbWFC_3d_all.tif','tiff','resolution',300)
 
 
 
 
 
 
 
 