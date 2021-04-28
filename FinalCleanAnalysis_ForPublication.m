clear
load CB1analysis cb1
load NMDAanalysis nmA
%This is to stract the data from the original structures. It makes it
%easier to analyze using ANOVA and other tests.
allV=[];
allG=[];
normalization='glnorm';%'l2l' %glnorm
for receptor='nc'
    if strcmp(receptor,'n')
        thisReceptor=nmA;
    else
        thisReceptor=cb1;
    end
    for layer={'ML','PL','GL'}
        for genetic={'WT','BT'}
            for env='CE'
                for gender='MF'
                    for animal=1:length(thisReceptor.(genetic{1}).(gender).(env))%aimal number
                        thisV=givemeV(thisReceptor.(genetic{1}).(gender).(env){animal},layer{1},normalization);
                        thisLabel=[genetic{1}(1) gender env num2str(animal) receptor lower(layer{1}(1))];
                        allV=[allV; thisV.values'];
                        allG=[allG; repmat(thisLabel,[3 1])];
                        
                    end
                end
            end
        end
    end
end
%% I just reorganize the data because of previous ways we analyzed them
avgG=[];
avgV=[];
newallV=[];newallG=[];avgSG=[];avgSV=[];
for d='CE'%standard or enriched environment
    for c='WB'%c57bl6 or BTBR
        for b='nc'%nmda or cb1
            for a='MF'%male or female
                thisBool=(allG(:,5)==b).*(allG(:,2)==a).*(allG(:,1)==c).*(allG(:,3)==d);
                thisG=allG(logical(thisBool),:);
                thisV=allV(logical(thisBool));
              
                %thisVml=mean(reshape(thisV,[9 3]),2);
                %thisGml=thisG(1:9,1:5);
%                 avgG=[avgG; thisGml];
%                 avgV=[avgV; thisVml];
                %the samples without averaging the 3 slices per animal.
                newallG=[newallG; thisG];
                newallV=[newallV; allV(logical(thisBool),:)];
                
%                 %average the 3 samples per animal
                 avg3=mean(reshape(thisV,[3,length(thisV)/3]))';
                 avgSV=[avgSV; avg3];
                 avgSG=[avgSG; thisG(1:3:end,:)];
              
            end
        end
    end
end


%% 
% Compare layers
% we are using newallG and newallV which have all the images (3 images per
% moouse)
%molecular
%This section was not used for the paper, but you can run it to see the
%numbers
figure(2)
thisBool=(newallG(:,6)=='m');
thisG=newallG(logical(thisBool),:);
thisV=newallV(logical(thisBool));
[sexTbML,groupstbML,gposthoctbML,figML,stdV]=  compareGroups(thisV,thisG,'anovan','Bonferroni');
writetable(sexTbML,'KW_genderML.xlsx')
writetable(groupstbML,'KW_grupsML.xlsx')
writetable(gposthoctbML,'PostHocGroupsML.xlsx')
figML.WindowStyle='normal';
set(gcf,'renderer','Painters','Units','inches')
exportgraphics(gcf,'MLAllBars.eps')
figML.WindowStyle='docked';


%Purkinje
figure(3)
thisBool=(newallG(:,6)=='p');
thisG=newallG(logical(thisBool),:);
thisV=newallV(logical(thisBool));
[sexTbPL,groupstbPL,gposthoctbPL,figPL]=  compareGroups(thisV,thisG,'anovan','Bonferroni');
writetable(sexTbPL,'KW_genderPL.xlsx')
writetable(groupstbPL,'KW_grupsPL.xlsx')
writetable(gposthoctbPL,'PostHocGroupsPL.xlsx')
set(gcf,'renderer','Painters','Units','inches')
exportgraphics(gcf,'PLAllBars.eps')

%Granular
figure(4)
thisBool=(newallG(:,6)=='g');
thisG=newallG(logical(thisBool),:);
thisV=newallV(logical(thisBool));
[sexTbGL,groupstbGL,gposthoctbGL,figGL]=  compareGroups(thisV,thisG,'anovan','Bonferroni');
writetable(sexTbGL,'KW_genderGL.xlsx')
writetable(groupstbGL,'KW_grupsGL.xlsx')
writetable(gposthoctbGL,'PostHocGroupsGL.xlsx')
set(gcf,'renderer','Painters','Units','inches')
exportgraphics(gcf,'GLAllBars.eps')

%% compare when averaging the slices from the same animal
%compare sex, gender, layer, environment
%This produces Figure 2

thisBool=(avgSG(:,5)=='c');%all cb1
thisG=avgSG(logical(thisBool),:);
thisV=avgSV(logical(thisBool));
[~,allCB1Tb,allCB1Tbst]=anovan(thisV,{thisG(:,1) thisG(:,2) thisG(:,3) thisG(:,6)},'model','linear','varnames',{'Gen','Sex','Env','Layer'},'Display','off');
dummyTb=array2table(allCB1Tb(2:end,1:end),'VariableNames',allCB1Tb(1,1:end));
eta2=[];
c=1;
for a=1:6
    eta2(c)=dummyTb{a,2}{1}./dummyTb{end,2}{1};
    c=c+1;
end    
dummyTb.eta2=eta2';
writetable(dummyTb,'MANOVA_allCB.xlsx')


[anovaCM,posthocCM]=comparegroupslayer(avgSV, avgSG,'c','M');
[anovaCF,posthocCF]=comparegroupslayer(avgSV, avgSG,'c','F');
     
writetable(anovaCM,'MANOVA_maleCB.xlsx')
writetable(anovaCF,'MANOVA_femCB.xlsx')


thisBool=(avgSG(:,5)=='n');%all nmda
thisG=avgSG(logical(thisBool),:);
thisV=avgSV(logical(thisBool));
[dummy,allNMTb,allNMTbst]=anovan(thisV,{thisG(:,1) thisG(:,2) thisG(:,3) thisG(:,6)},'model','linear','varnames',{'Gen','Sex','Env','Layer'},'Display','off');
dummyTb=array2table(allNMTb(2:end,:),'VariableNames',allNMTb(1,:));
eta2=[];
c=1;
for a=1:6
    eta2(c)=dummyTb{a,2}{1}./dummyTb{end,2}{1};
    c=c+1;
end    
dummyTb.eta2=eta2';
writetable(dummyTb,'MANOVA_allNM.xlsx')

[anovaNM,posthocNM]=comparegroupslayer(avgSV, avgSG,'n','M');
[anovaNF,posthocNF]=comparegroupslayer(avgSV, avgSG,'n','F');
     
writetable(anovaNM,'MANOVA_maleNM.xlsx')
writetable(anovaNF,'MANOVA_femNM.xlsx')

%average the values of all the layers. 
[anovaNMavl,posthocNMavl,groupnames,nMval]=comparegroupsAvlayer(avgSV, avgSG,'n','M');
[anovaNFavl,posthocNFavl,~,nFval]=comparegroupsAvlayer(avgSV, avgSG,'n','F');
[anovaCMavl,posthocCMavl,groupnames,cMval]=comparegroupsAvlayer(avgSV, avgSG,'c','M');
[anovaCFavl,posthocCFavl,~,cFval]=comparegroupsAvlayer(avgSV, avgSG,'c','F');

figure(5)
clf
subplot 221
plotexpbars(nMval,groupnames)

subplot 222
plotexpbars(nFval,groupnames)

subplot 223
plotexpbars(cMval,groupnames)

subplot 224
plotexpbars(cFval,groupnames)
figh=gcf;
figh.WindowStyle='normal';
set(figh,'renderer','Painters','Units','inches')
set(figh,'WindowStyle','normal')
set(figh, 'PaperUnits', 'inches','Renderer','painters');
set(figh,'Units','inches','Position',[1 1 4 4])
exportgraphics(gcf,'FinalAnalysis.eps')
figh.WindowStyle='docked';

writetable([anovaNMavl; anovaNFavl;anovaCMavl;anovaCFavl],'MANOVA_avLayerAll.xlsx')
writetable([posthocNMavl],'posthoc_avLayerAll.xlsx','Sheet',1);
writetable([posthocNFavl],'posthoc_avLayerAll.xlsx','Sheet',2,'WriteMode','append');
writetable([posthocCMavl],'posthoc_avLayerAll.xlsx','Sheet',3,'WriteMode','append');
writetable([posthocCFavl],'posthoc_avLayerAll.xlsx','Sheet',4,'WriteMode','append');

%% Not used for the paper
%just for curiosity compare per layer
% This comparision suffers from having too few data points, but the trends
% remain.

[anovaNMm,posthocNMm,groupnames,nMvalm]=comparegroupsOnelayer(avgSV, avgSG,'n','M','m');
[anovaNMp,posthocNMp,groupnames,nMvalp]=comparegroupsOnelayer(avgSV, avgSG,'n','M','p');
[anovaNMg,posthocNMg,groupnames,nMvalg]=comparegroupsOnelayer(avgSV, avgSG,'n','M','g');

[anovaNFm,posthocNFm,groupnames,nFvalm]=comparegroupsOnelayer(avgSV, avgSG,'n','F','m');
[anovaNFp,posthocNFp,groupnames,nFvalp]=comparegroupsOnelayer(avgSV, avgSG,'n','F','p');
[anovaNFg,posthocNFg,groupnames,nFvalg]=comparegroupsOnelayer(avgSV, avgSG,'n','F','g');

[anovaCMm,posthocCMm,groupnames,cMvalm]=comparegroupsOnelayer(avgSV, avgSG,'c','M','m');
[anovaCMp,posthocCMp,groupnames,cMvalp]=comparegroupsOnelayer(avgSV, avgSG,'c','M','p');
[anovaCMg,posthocCMg,groupnames,cMvalg]=comparegroupsOnelayer(avgSV, avgSG,'c','M','g');

[anovaCFm,posthocCFm,groupnames,cFvalm]=comparegroupsOnelayer(avgSV, avgSG,'c','F','m');
[anovaCFp,posthocCFp,groupnames,cFvalp]=comparegroupsOnelayer(avgSV, avgSG,'c','F','p');
[anovaCFg,posthocCFg,groupnames,cFvalg]=comparegroupsOnelayer(avgSV, avgSG,'c','F','g');

figure(6)
clf
subplot 231
plotexpbars(nMvalm,groupnames)
subplot 232
plotexpbars(nMvalp,groupnames)
subplot 233
plotexpbars(nMvalg,groupnames)
subplot 234
plotexpbars(nFvalm,groupnames)
subplot 235
plotexpbars(nFvalp,groupnames)
subplot 236
plotexpbars(nFvalg,groupnames)

figure(7)
clf
subplot 231
plotexpbars(cMvalm,groupnames)
subplot 232
plotexpbars(cMvalp,groupnames)
subplot 233
plotexpbars(cMvalg,groupnames)
subplot 234
plotexpbars(cFvalm,groupnames)
subplot 235
plotexpbars(cFvalp,groupnames)
subplot 236
plotexpbars(cFvalg,groupnames)



%% lillifors
% this checks for the normal distribution of all experimental groups.
% At the end I calculate the number of sets that are not normal. 
layer={'m','p','g'};
rec={'n','c'};
ee={'C','E'};
sex={'M','F'};
gen={'W','B'};

c1=1;
for a=layer
    for b=rec
        for c=ee
            for d=sex
                for e=gen
                    myG(c1,:)=[e{1} d{1} c{1} b{1} a{1}];
                    boolV=(allG(:,1)==e{1}).*(allG(:,2)==d{1}).*(allG(:,3)==c{1}).*(allG(:,5)==b{1}).*(allG(:,6)==a{1});
                    thisV=allV(logical(boolV));
                    thisG=allG(logical(boolV),:);
                    [h,p,ks,critval]=lillietest(thisV(1:end));
                    llt(c1)=h;
                    llp(c1)=round(100*p)/100;
                    ksV(c1)=ks;
                    cvV(c1)=p;
                    disp(['hip ', num2str(h), ' p ',num2str(p)  ' ',a{1}, b{1}, c{1}, d{1}, e{1}])
                    c1=c1+1;
                end
            end
        end
    end
end
nnz(llt)

%% 
%% Compare slices from the same animal
%for simplicity I am using the values for the granule cell layers in NMDAR
%and CB1 animals, both male and females. 
boolV=(allG(:,6)=='g').*(allG(:,5)=='n');
nmV=allV(logical(boolV));
nmG=allG(logical(boolV),1:4);
boolV=(allG(:,6)=='g').*(allG(:,5)=='c');
cbV=allV(logical(boolV));
cbG=allG(logical(boolV),1:4);
[fV,fG]=compareSameAn(nmV,nmG,cbV,cbG);
groups={'WMC','WME','BMC','BME','WFC','WFE','BFC','BFE'};
groupsLabel={'B6MC','B6ME','BTMC','BTME','B6FC','B6FE','BTFC','BTFE'};
groupsLabel={'C57s','C57e','BTs','BTe','C57s','C57e','BTs','BTe'};


figure(3)
clf
for a=1:length(groups)
    bool=(sum((fG(:,1:3)==repmat(groups{a},[length(fG) 1])),2)==3);
    boolV=fV(logical(bool),:);
    for b=1:3
     
      hold on    
      if groups{a}(2)=='M'
          h=scatter(boolV(b,2),boolV(b,1),350,'d');
          set(h,'Marker','o','SizeData',550,'MarkerFaceColor',[0.7 0.7 0.7],...
              'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceAlpha',0.7);
          h=text(boolV(b,2)-0.25,boolV(b,1),groupsLabel{a});
      else
          h=scatter(boolV(b,2),boolV(b,1),350,'d');
          set(h,'Marker','sq','SizeData',750,'MarkerEdgeColor','k')
          h=text(boolV(b,2)-0.32,boolV(b,1),groupsLabel{a});
      end
      
      set(h,'FontSize',8)
    end
end
axis([0 7 0 15])
axis square
ylabel('NMDAR1/DAPI')
xlabel('CB1R/DAPI')

figh=gcf;
figh.WindowStyle='normal';
set(figh,'renderer','Painters','Units','inches')
set(figh,'WindowStyle','normal')
set(figh, 'PaperUnits', 'inches','Renderer','painters');
set(figh,'Units','inches','Position',[1 1 4 4])
exportgraphics(gcf,'ScatterPlot.eps')
figh.WindowStyle='docked';



boolV=(fG(:,2)=='M');
fV2=fV(logical(boolV),:);
[r,p]=corrcoef(fV2(:,1),fV2(:,2))

boolV=(fG(:,2)=='F');
fV2=fV(logical(boolV),:);
[r,p]=corrcoef(fV2(:,1),fV2(:,2))

boolV=(fG(:,1)=='B').*(fG(:,2)=='M');
fV2=fV(logical(boolV),:);
[r,p]=corrcoef(fV2(:,1),fV2(:,2))

boolV=(fG(:,1)=='B').*(fG(:,2)=='F');
fV2=fV(logical(boolV),:);
[r,p]=corrcoef(fV2(:,1),fV2(:,2))

%% Images showing the ROI example in the Supplementary Figure S1
clear
nWMC=load('cbWCM3.mat','cb');
nWME=load('cbWEM3.mat','cb');
nWFC=load('cbWCF2.mat','cb');
nWFE=load('cbWEF.mat','cb');
nBMC=load('cbBCM.mat','cb');
nBME=load('cbBEM3.mat','cb');
nBFC=load('cbBCF.mat','cb');
nBFE=load('cbBEF2.mat','cb');

thisSt=nWMC.cb.WT;%
an2a=1;
im2a=3;
rcolor=2;
thisR=thisSt.R.M.C{an2a}(im2a);
thisD=thisSt.D.M.C{an2a}(im2a);
bcxy=thisR.Data.bcxy;
thisroiSt=thisR.ROI;
thisRim=thisR.Data.im;
thisDim=thisD.Data.im;
rgbM=zeros([size(thisRim) 3]);
rgbM(:,:,rcolor)=thisRim./max(thisRim(:));
rgbM(:,:,3)=thisDim./max(thisDim(:));
mlroi=round(thisroiSt.Layer(1).ROIxy);
plroi=thisroiSt.Layer(2).ROIxy;
glroi=thisroiSt.Layer(3).ROIxy;


clf
imagesc(rgbM)
axis equal
axis off
hold on
polml=polyshape(mlroi(:,1),mlroi(:,2));
plot(polml)
polgl=polyshape(glroi(:,1),glroi(:,2));
plot(polgl)
pcroi=[-25 -25;-25 25; 25 25; 25 -25]; %from newROI.m I know, this is terrible programming.
for a=1:size(plroi,1)
    polpl=polyshape(plroi(a,1)+pcroi(:,1),plroi(a,2)+pcroi(:,2));
    plot(polpl)
end
figML=gcf;
figML.WindowStyle='normal';
set(gcf,'renderer','Painters','Units','inches')
exportgraphics(gca,'ExampleROI.tif')
figML.WindowStyle='docked';
dx=0.388;


cb{1}=nWMC.cb.WT.R.M.C{1}(3:5);
cb{2}=nWMC.cb.WT.R.M.C{2}(2:4);
cb{3}=nWMC.cb.WT.R.M.C{3}(1:3);
[glV(:,:,1),plV(:,:,1),mlV(:,:,1),pcV(:,:,1)]=getmeROIstats(cb);
cb{1}=nWME.cb.WT.R.M.E{1}([4 5 7]); 
cb{2}=nWME.cb.WT.R.M.E{2}(3:5);
cb{3}=nWME.cb.WT.R.M.E{3}(1:3);
[glV(:,:,2),plV(:,:,2),mlV(:,:,2),pcV(:,:,2)]=getmeROIstats(cb);
cb{1}=nWFC.cb.WT.R.F.C{2}([2 3 4]); 
cb{2}=nWFC.cb.WT.R.F.C{3}([3 7:8]);
cb{3}=nWFC.cb.WT.R.F.C{4}(3:5);
[glV(:,:,3),plV(:,:,3),mlV(:,:,3),pcV(:,:,3)]=getmeROIstats(cb);
cb{1}=nWFE.cb.WT.R.F.E{2}(7:9);
cb{2}=nWFE.cb.WT.R.F.E{3}(1:3);
cb{3}=nWFE.cb.WT.R.F.E{4}(6:8);
[glV(:,:,4),plV(:,:,4),mlV(:,:,4),pcV(:,:,4)]=getmeROIstats(cb);
cb{1}=nBMC.cb.BT.R.M.C{1}(2:4); 
cb{2}=nBMC.cb.BT.R.M.C{2}(7:9);
cb{3}=nBMC.cb.BT.R.M.C{3}([4 5 8]);
[glV(:,:,5),plV(:,:,5),mlV(:,:,5),pcV(:,:,5)]=getmeROIstats(cb);
cb{1}=nBME.cb.BT.R.M.E{2}(2:4);%red0
cb{2}=nBME.cb.BT.R.M.E{3}(1:3);
cb{3}=nBME.cb.BT.R.M.E{4}(3:5);
[glV(:,:,6),plV(:,:,6),mlV(:,:,6),pcV(:,:,6)]=getmeROIstats(cb);  
cb{1}=nBFC.cb.BT.R.F.C{1}(5:7);
cb{2}=nBFC.cb.BT.R.F.C{2}([3 4 8]);
cb{3}=nBFC.cb.BT.R.F.C{3}(1:3);
[glV(:,:,7),plV(:,:,7),mlV(:,:,7),pcV(:,:,7)]=getmeROIstats(cb);
 cb{1}=nBFE.cb.BT.R.F.E{2}([1 3 5]);
 cb{2}=nBFE.cb.BT.R.F.E{3}(3:5);
  cb{3}=nBFE.cb.BT.R.F.E{4}(1:3);
[glV(:,:,8),plV(:,:,8),mlV(:,:,8),pcV(:,:,8)]=getmeROIstats(cb);


%%
% %% Old stuff not used for the paper anymore
% %compare layers
% %molecular
% figure(5)
% thisBool=(avgSG(:,6)=='m');
% thisG=avgSG(logical(thisBool),:);
% thisV=avgSV(logical(thisBool));
% [sexTbML,groupstbML,gposthoctbML,figML]=  compareGroups(thisV,thisG,'anovan','Bonferroni');
% writetable(sexTbML,'KW_genderMLavgSamples.xlsx')
% writetable(groupstbML,'KW_grupsMLavgSamples.xlsx')
% writetable(gposthoctbML,'PostHocGroupsMLavgSamples.xlsx')
% figML.WindowStyle='normal';
% set(gcf,'renderer','Painters','Units','inches')
% exportgraphics(gcf,'MLAllBarsAvgSamples.eps')
% figML.WindowStyle='docked';
% 
% figure(6)
% thisBool=(avgSG(:,6)=='p');
% thisG=avgSG(logical(thisBool),:);
% thisV=avgSV(logical(thisBool));
% [sexTbML,groupstbML,gposthoctbML,figML]=  compareGroups(thisV,thisG,'anovan','');
% writetable(sexTbML,'KW_genderPLavgSamples.xlsx')
% writetable(groupstbML,'KW_grupsPLavgSamples.xlsx')
% writetable(gposthoctbML,'PostHocGroupsPLavgSamples.xlsx')
% figML.WindowStyle='normal';
% set(gcf,'renderer','Painters','Units','inches')
% exportgraphics(gcf,'PLAllBarsAvgSamples.eps')
% figML.WindowStyle='docked';
% 
% figure(7)
% thisBool=(avgSG(:,6)=='g');
% thisG=avgSG(logical(thisBool),:);
% thisV=avgSV(logical(thisBool));
% [sexTbML,groupstbML,gposthoctbML,figML]=  compareGroups(thisV,thisG,'anovan','');
% writetable(sexTbML,'KW_genderGLavgSamples.xlsx')
% writetable(groupstbML,'KW_grupsGLavgSamples.xlsx')
% writetable(gposthoctbML,'PostHocGroupsGLavgSamples.xlsx')
% figML.WindowStyle='normal';
% set(gcf,'renderer','Painters','Units','inches')
% exportgraphics(gcf,'GLAllBarsAvgSamples.eps')
% figML.WindowStyle='docked';
% 
