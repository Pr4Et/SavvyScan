function normpic2=balance_pic(normpic,Nshades)    
%Written by Shahar Seifer, Weizmann Institute of Science, 2020
%H=histogram(normpic,'BinMethod','auto');
[BinValues,BinEdges]=histcounts(normpic);
NumBins=length(BinValues);    
sumH=sum(BinValues);
    temp=0;
    for n=1:NumBins-1
        temp=temp+BinValues(n);
        if temp>0.005*sumH
            lowedge=BinEdges(n);
        break;
        end
    end
    temp=0;
    for n2=NumBins:-1:2
        temp=temp+BinValues(n2);
        if temp>0.005*sumH
            highedge=BinEdges(n2);
        break;
        end
    end
    normpic(normpic>highedge)=highedge; %remove white dots
    normpic(normpic<lowedge)=lowedge; %remove black dots
    normpic2=double(normpic-lowedge)*double(Nshades)/(highedge-lowedge);
return 
    