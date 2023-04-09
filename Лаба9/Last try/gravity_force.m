function fe1 = gravity_force(ige1m, ige2m, xc, yc, assoc, els, fe)
    d1=2300;
    for i=1:size(ige1m,2)
       tmp = ige1m(i);
       xi = xc(assoc(tmp,1));
       yi = yc(assoc(tmp,1));
       
       xj = xc(assoc(tmp,2));
       yj = yc(assoc(tmp,2));
       
       xk = xc(assoc(tmp,3));
       yk = yc(assoc(tmp,3));
       
       S = get_square(xi, yi, xj, yj, xk ,yk);
       
       it = assoc(tmp,1);
       jt = assoc(tmp,2);
       kt = assoc(tmp,3);
       if(yi~=0)
        fe(2*it) = fe(2*it) - d1*9.8*S/3;
       end
       if(yj~=0)
        fe(2*jt) = fe(2*jt) - d1*9.8*S/3;
       end
       if(yk~=0)
       fe(2*kt) = fe(2*kt) - d1*9.8*S/3;
       end
    end
    
    for i=1:size(ige2m,2)
       tmp = ige2m(i);
       xi = xc(assoc(tmp,1));
       yi = yc(assoc(tmp,1));
       
       xj = xc(assoc(tmp,2));
       yj = yc(assoc(tmp,2));
       
       xk = xc(assoc(tmp,3));
       yk = yc(assoc(tmp,3));
       
       S = get_square(xi, yi, xj, yj, xk ,yk);
       
       it = assoc(tmp,1);
       jt = assoc(tmp,2);
       kt = assoc(tmp,3);
       if(yi~=0)
        fe(2*it) = fe(2*it) - d1*9.8*S/3;
       end
       if(yj~=0)
        fe(2*jt) = fe(2*jt) - d1*9.8*S/3;
       end
       if(yk~=0)
       fe(2*kt) = fe(2*kt) - d1*9.8*S/3;
       end
    end
    fe1=fe;
end

