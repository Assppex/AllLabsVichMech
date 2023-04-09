function force = NPUf(xnds,ynds,nds,ige1,base)

kolnds = size(ige1,2);
underload = [];
force = zeros(2*nds,1);
for i=1:kolnds
   curnd = ige1(i);
   
   x = xnds(curnd);
   y = ynds(curnd);

   if(abs(y-275)<=10^(-4) && x <= 424.7 && x>= 412.3)
      underload(end+1) = curnd;
   end
    
end
forceNPUunsorted = underload;
% сортируем узлы так, чтобы они стояли слевуа направо
for j=1:size(underload,2)-1
   for k=1:size(underload,2)-j
       x = xnds(underload(k));
       xnxt = xnds(underload(k+1));
      if(x > xnxt)
          tmp = underload(k+1);
          underload(k+1) = underload(k);
          underload(k) = tmp;
          
      end
   end
end
size(underload,2)
for k=1:size(underload,2)-1
    x = xnds(underload(k));
    xnxt = xnds(underload(k+1));
    
    y = ynds(underload(k));
    ynxt = ynds(underload(k+1));
    
    length = sqrt((xnxt-x)^2 + (ynxt-y)^2);
   
   force(2*underload(k)) = force(2*underload(k)) - 500*9.8*15*length;
   force(2*underload(k+1)) = force(2*underload(k+1)) - 500*9.8*15*length;
   
end



kolnds = size(ige1,2);
underload = [];
for i=1:kolnds
   curnd = ige1(i);
   
   x = xnds(curnd);
   y = ynds(curnd);

   if(abs(y-2/5*x-105.2)<=10^(-4) && x <=462  && x >= 424.5 )
      underload(end+1) = curnd;
   end
    
end
forceNPUunsorted = underload;
% сортируем узлы так, чтобы они стояли слевуа направо
for j=1:size(underload,2)-1
   for k=1:size(underload,2)-j
       x = xnds(underload(k));
       xnxt = xnds(underload(k+1));
      if(x > xnxt)
          tmp = underload(k+1);
          underload(k+1) = underload(k);
          underload(k) = tmp;
          
      end
   end
end

for k=1:size(underload,2)-1
    x = xnds(underload(k));
    xnxt = xnds(underload(k+1));
    
    y = ynds(underload(k));
    ynxt = ynds(underload(k+1));
    
    length = sqrt((xnxt-x)^2 + (ynxt-y)^2);
   
   force(2*underload(k)) = force(2*underload(k)) - 500*9.8*(290-y)*length*cos(atan(2/5));
   force(2*underload(k)+1) = force(2*underload(k)+1) + 500*9.8*(290-y)*length*sin(atan(2/5));
   
   force(2*underload(k+1)) = force(2*underload(k+1)) - 500*9.8*(290-ynxt)*length*cos(atan(2/5));
   force(2*underload(k+1)+1) = force(2*underload(k+1)+1) + 500*9.8*(290-ynxt)*length*sin(atan(2/5));
   
   
end







kolnds = size(base,2);
underload = [];
for i=1:kolnds
   curnd = base(i);
   
   x = xnds(curnd);
   y = ynds(curnd);

   if(y==0 && x <= 0)
      underload(end+1) = curnd;
   end
    
end
forceNPUunsorted = underload;
% сортируем узлы так, чтобы они стояли слевуа направо
for j=1:size(underload,2)-1
   for k=1:size(underload,2)-j
       x = xnds(underload(k));
       xnxt = xnds(underload(k+1));
      if(x > xnxt)
          tmp = underload(k+1);
          underload(k+1) = underload(k);
          underload(k) = tmp;
          
      end
   end
end
size(underload,2)
for k=1:size(underload,2)-1
    x = xnds(underload(k));
    xnxt = xnds(underload(k+1));
    
    y = ynds(underload(k));
    ynxt = ynds(underload(k+1));
    
    length = sqrt((xnxt-x)^2 + (ynxt-y)^2);
   
   force(2*underload(k)) = force(2*underload(k)) - 500*9.8*290*length;
   force(2*underload(k+1)) = force(2*underload(k+1)) - 500*9.8*290*length;
   
end


kolnds = size(ige1,2);
underload = [];
for i=1:kolnds
   curnd = ige1(i);
   
   x = xnds(curnd);
   y = ynds(curnd);

   if(abs(y-2/3*x)<=10^(-4) && x >= 0 && x <= 412.5)
      underload(end+1) = curnd;
   end
    
end

% сортируем узлы так, чтобы они стояли слевуа направо
for j=1:size(underload,2)-1
   for k=1:size(underload,2)-j
       x = xnds(underload(k));
       xnxt = xnds(underload(k+1));
      if(x > xnxt)
          tmp = underload(k+1);
          underload(k+1) = underload(k);
          underload(k) = tmp;
          
      end
   end
end
size(underload,2)
for k=1:size(underload,2)-1
    x = xnds(underload(k));
    xnxt = xnds(underload(k+1));
    
    y = ynds(underload(k));
    ynxt = ynds(underload(k+1));
    
    length = sqrt((xnxt-x)^2 + (ynxt-y)^2);
    alpha = atan(2/3);
   
   force(2*underload(k)) = force(2*underload(k)) - 500*9.8*(290-y)*length*cos(alpha);
   force(2*underload(k)+1) = force(2*underload(k)+1) + 500*9.8*(290-y)*length*sin(alpha);
   
   force(2*underload(k+1)) = force(2*underload(k+1)) - 500*9.8*(290-ynxt)*length*cos(alpha);
   force(2*underload(k+1)+1) = force(2*underload(k+1)+1) + 500*9.8*(290-ynxt)*length*sin(alpha);
   
end


end

