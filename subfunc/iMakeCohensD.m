function cohenD = iMakeCohensD(data1, data2)

if size(data1,1) == 1 || size(data1,2) == 1 %if it's a vector
    pooledFun = @(s1,s2,n1,n2) sqrt(((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2));
    
    s1=nanstd(data1);
    s2=nanstd(data2);
    n1 =sum(~isnan(data1));
    n2 =sum(~isnan(data2));
    
    xbar1 = nanmean(data1);
    xbar2 = nanmean(data2);
    
    cohenD=(xbar2-xbar1)/pooledFun(s1,s2,n1,n2);
    
else
    pooledFun = @(s1,s2,n1,n2) sqrt(((n1-1).*s1.^2+(n2-1).*s2.^2)./(n1+n2-2));
    
    s1=nanstd(data1);
    s2=nanstd(data2);
    n1 =sum(~isnan(data1));
    n2 =sum(~isnan(data2));
    
    xbar1 = nanmean(data1);
    xbar2 = nanmean(data2);
    
    cohenD=(xbar2-xbar1)./pooledFun(s1,s2,n1,n2);
end