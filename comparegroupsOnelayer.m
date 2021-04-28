function [dummyTb,gposthoctb,gname,resV]=comparegroupsOnelayer(avgSV, avgSG,receptor,sex,layer)


thisBool=((avgSG(:,5)==receptor).*(avgSG(:,2)==sex).*(avgSG(:,6)==layer));
thisG=avgSG(logical(thisBool),:);
thisV=avgSV(logical(thisBool));
%average of all the layers
thisV1=thisV;
thisG1=thisG;

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
resV=reshape(thisV1,[3, 4]);
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
