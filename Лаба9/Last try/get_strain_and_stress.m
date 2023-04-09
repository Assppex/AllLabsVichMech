function [strain, stress] = get_strain_and_stress(ux,uy, els, xnds, ynds, assoc, ige1m, ige2m,basem,E1,v1,E2,v2,E3,v3)
    straintmp = zeros(3, els);
    stresstmp = zeros(3, els);
    for k=1:size(ige1m,2)
       i = ige1m(k);
       it = assoc(i,1); 
       jt = assoc(i,2);
       kt = assoc(i,3);
       
       xi = xnds(it);
       yi = ynds(it);
       
       xj = xnds(jt);
       yj = ynds(jt);
       
       xk = xnds(kt);
       yk = ynds(kt);
       
       
       
       uxi = ux(it);
       uyi = uy(it);
       
       uxj = ux(jt);
       uyj = uy(jt);
       
       uxk = ux(kt);
       uyk = uy(kt);
       
       Displ = [uxi; uyi;uxj;uyj;uxk;uyk];
       
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
    D = E1*(1-v1)/((1+v1)*(1-2*v1))*[[1,v1/(1-v1),0];[v1/(1-v1),1,0];[0,0,(1-2*v1)/(2*(1-v1))]];
    straintmp(:,i) = B * Displ;
    stresstmp(:,i) = D*straintmp(:,i);
    end
    
    for k=1:size(ige2m,2)
       i = ige2m(k);
       it = assoc(i,1); 
       jt = assoc(i,2);
       kt = assoc(i,3);
       
       xi = xnds(it);
       yi = ynds(it);
       
       xj = xnds(jt);
       yj = ynds(jt);
       
       xk = xnds(kt);
       yk = ynds(kt);
       
       
       
       uxi = ux(it);
       uyi = uy(it);
       
       uxj = ux(jt);
       uyj = uy(jt);
       
       uxk = ux(kt);
       uyk = uy(kt);
       
       Displ = [uxi; uyi;uxj;uyj;uxk;uyk];
       
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
    D=E2*(1-v2)/((1+v2)*(1-2*v2))*[[1,v2/(1-v2),0];[v2/(1-v2),1,0];[0,0,(1-2*v2)/(2*(1-v2))]];
    straintmp(:,i) = B * Displ;
    stresstmp(:,i) = D*straintmp(:,i);
    end
    
    
    for k=1:size(basem,2)
       i = basem(k);
       it = assoc(i,1); 
       jt = assoc(i,2);
       kt = assoc(i,3);
       
       xi = xnds(it);
       yi = ynds(it);
       
       xj = xnds(jt);
       yj = ynds(jt);
       
       xk = xnds(kt);
       yk = ynds(kt);
       
       
       
       uxi = ux(it);
       uyi = uy(it);
       
       uxj = ux(jt);
       uyj = uy(jt);
       
       uxk = ux(kt);
       uyk = uy(kt);
       
       Displ = [uxi; uyi;uxj;uyj;uxk;uyk];
       
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
    D=E3*(1-v3)/((1+v3)*(1-2*v3))*[[1,v3/(1-v3),0];[v3/(1-v3),1,0];[0,0,(1-2*v3)/(2*(1-v3))]];
    straintmp(:,i) = B * Displ;
    stresstmp(:,i) = D*straintmp(:,i);
    end
    strain = straintmp(1:2,:)';
    stress = stresstmp(1:2,:)';
end

