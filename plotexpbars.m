function plotexpbars(cMval,groupnames)

ceMstd=std(cMval)/sqrt(size(cMval,1));
groupnames={'C57 S','BT S','C57 E','BT E'};

g1=bar([1 2 3 4],mean(cMval(:,[1 3 2 4])));
hold on
g2=plot([1 2 3 4],cMval(:,[1 3 2 4]),'ko');
g5=errorbar([1 2 3 4],mean(cMval(:,[ 1 3 2 4])),ceMstd([1 3 2 4]),'k.');
g1.FaceAlpha=0;g1.XData;
g5.LineWidth=1;g5.CapSize=10;g2(1).MarkerSize=3;g2(2).MarkerSize=3;g2(3).MarkerSize=3;
box off
pg1=g1.Parent;
pg1.XTickLabel=groupnames([1 3 2 4]);
pg1.YLim=[0 17];
pg1.XLim=[0.5 4.5];
