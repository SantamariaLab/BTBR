function [dummyTb,gposthoctb,gname,resV]=comparegroupsAvlayer(avgSV, avgSG,receptor,sex)


thisBool=((avgSG(:,5)==receptor).*(avgSG(:,2)==sex));%all cb1
thisG=avgSG(logical(thisBool),:);
thisV=avgSV(logical(thisBool));
%average of all the layers
thisV1=[mean(reshape(thisV(1:3:end,:),[3 4])) mean(reshape(thisV(2:3:end,:),[3 4])) mean(reshape(thisV(3:3:end,:),[3 4]))];
thisG1=[thisG(1:9:end,1:5);thisG(2:9:end,1:5);thisG(3:9:end,1:5)];

[~,maleCB1Tb,maleCB1Tbst]=anovan(thisV,{thisG(:,1) thisG(:,3) },'model','linear','varnames',{'Gen','Env'},'Display','off');
dummyTb=array2table(maleCB1Tb(2:end,:),'VariableNames',maleCB1Tb(1,:));
eta2=[];
c=1;
for a=1:4
    eta2(c)=dummyTb{a,2}{1}./dummyTb{end,2}{1};
    c=c+1;
end
dummyTb.eta2=eta2';
%writetable(dummyTb,'MANOVA_maleCB.xlsx')
[multi1,msd1,~,gname]=multcompare(maleCB1Tbst,'Dimension',[1 2],'Display','on');
resV=[thisV1(1:4:end)' thisV1(2:4:end)' thisV1(3:4:end)' thisV1(4:4:end)'];
meanV=mean(resV);
stdV=std(resV);
cohend=[(meanV(1)-meanV(2))./sqrt((stdV(1).^2+stdV(2).^2)/2);...
    (meanV(1)-meanV(3))./sqrt((stdV(1).^2+stdV(3).^2)/2);...
    (meanV(1)-meanV(4))./sqrt((stdV(1).^2+stdV(4).^2)/2);...
    (meanV(2)-meanV(3))./sqrt((stdV(2).^2+stdV(3).^2)/2);...
    (meanV(2)-meanV(4))./sqrt((stdV(2).^2+stdV(4).^2)/2);...
    (meanV(3)-meanV(4))./sqrt((stdV(3).^2+stdV(4).^2)/2)];
gpostArray=[gname(multi1(:,1)) gname(multi1(:,2)) num2cell(multi1(:,6)) num2cell(cohend(:))];
gposthoctb=array2table(gpostArray,...
    'VariableNames',[{'group1'},{'group2'},{[receptor sex]},{'cohen'}]);
