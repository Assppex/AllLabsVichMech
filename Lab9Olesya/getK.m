function getK = getK(ige1m,ige2m,basem, assoc, els, xnds, ynds, E1,v1,E2,v2,E3,v3,nds)
K = zeros(2*nds,2*nds);
for i=1:els
   for j=1:size(ige1m,2)
      if(i==ige1m(j))
         D = E1*(1-v1)/((1+v1)*(1-2*v1))*[[1,v1/(1-v1),0];[v1/(1-v1),1,0];[0,0,(1-2*v1)/(2*(1-v1))]];
         break;
      end
   end
   
   for j=1:size(ige2m,2)
      if(i==ige2m(j))
         D=E2*(1-v2)/((1+v2)*(1-2*v2))*[[1,v2/(1-v2),0];[v2/(1-v2),1,0];[0,0,(1-2*v2)/(2*(1-v2))]];
         break;
      end
   end
   
   for j=1:size(basem,2)
      if(i==basem(j))
         D=E3*(1-v3)/((1+v3)*(1-2*v3))*[[1,v3/(1-v3),0];[v3/(1-v3),1,0];[0,0,(1-2*v3)/(2*(1-v3))]];
         break;
      end
   end
   xi=xnds(assoc(i,1));
    yi=ynds(assoc(i,1));
    
    xj=xnds(assoc(i,2));
    yj=ynds(assoc(i,2));
    
    xk=xnds(assoc(i,3));
    yk=ynds(assoc(i,3));
   der = [[-1,0,1];[-1,1,0]];
    J = der*[[xi,yi];[xj,yj];[xk,yk]];
    Ji=inv(J);
%     считаем В
    B=zeros(3,6);
    tmp = zeros(2,1);
    tmp = Ji*der(:,1);
    B(1,1)=tmp(1);
    B(3,1)=tmp(2);
    B(2,2)=tmp(2);
    B(3,2)=tmp(1);
    
    tmp = Ji*der(:,2);
    B(1,3)=tmp(1);
    B(3,3)=tmp(2);
    B(2,4)=tmp(2);
    B(3,4)=tmp(1);
    
    tmp = Ji*der(:,3);
    B(1,5)=tmp(1);
    B(3,5)=tmp(2);
    
    B(2,6)=tmp(2);
    B(3,6)=tmp(1);
    
    S = get_square(xi,yi,xj,yj,xk,yk);
    
    ke=B'*D*B*S;
    
    A=zeros(6,2*nds);
    
    it=assoc(i,1);
    jt=assoc(i,2);
    kt=assoc(i,3);
    
    A(1,2*it-1)=1;
    A(2,2*it)=1;
    A(3,2*jt-1)=1;
    A(4,2*jt)=1;
    A(5,2*kt-1)=1;
    A(6,2*kt)=1;
    
    K=K+A'*ke*A;
    
end
getK=K;
end

