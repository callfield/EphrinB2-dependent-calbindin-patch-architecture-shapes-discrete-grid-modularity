function corr = nancorr2(a,b)
a = a - nanmean(a(:));
b = b - nanmean(b(:));
a(isnan(a))=0;
b(isnan(b))=0;
corr = sum(sum(a.*b))/sqrt(sum(sum(a.*a))*sum(sum(b.*b)));
end