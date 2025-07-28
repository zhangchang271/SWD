function p=pad(p0,nbc,isfs,pad_top)

if (isfs)
   p=[repmat(p0(:,1),1,nbc),p0,repmat(p0(:,end),1,nbc)];
   p=[repmat(p(1,:),pad_top,1);p;repmat(p(end,:),nbc,1)];
else    
   p=[repmat(p0(:,1),1,nbc),p0,repmat(p0(:,end),1,nbc)];
   p=[repmat(p(1,:),nbc,1);p;repmat(p(end,:),nbc,1)];
end 