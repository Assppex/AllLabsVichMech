function fe1 = gravity_force(ige1, ige2, xc, yc, assoc, els, fe)
        d1=2300;
%     for i=1:size(ige1m,2)
% 
%         a = ige1m(i);
%         xi = xc(assoc(a, 1));
%         yi = xc(assoc(a, 1));
% 
%         xj = xc(assoc(a, 2));
%         yj = xc(assoc(a, 2));
% 
%         xk = xc(assoc(a, 3));
%         yk = xc(assoc(a, 3));
%         
%         S = get_square(xi, yi, xj, yj, xk, yk);
%         
%         fe(2*a)=fe(2*a)-d1*9.8*S;
%     end
%     
%     for i=1:size(ige2,2)
%         a = ige2m(i);
%         xi = xc(assoc(a, 1));
%         yi = xc(assoc(a, 1));
% 
%         xj = xc(assoc(a, 2));
%         yj = xc(assoc(a, 2));
% 
%         xk = xc(assoc(a, 3));
%         yk = xc(assoc(a, 3));
%         
%         S = get_square(xi, yi, xj, yj, xk, yk);
%         S
%         fe(2*a)=fe(2*a)-d1*9.8*S;
%     end
    
       for i=1:size(ige1,2)
        tmp = [];
        a = ige1(i);
        for j=1:els
            it = assoc(j, 1);
            jt = assoc(j, 2);
            kt = assoc(j, 3);
            if(a == it || a==jt || a==kt)
                tmp = [tmp j];
            end
        end
        V = 0;
        for k=1:size(tmp,2)
            xi = xc(assoc(k, 1));
            yi = yc(assoc(k, 1));

            xj = xc(assoc(k, 2));
            yj = yc(assoc(k, 2));

            xk = xc(assoc(k, 3));
            yk = yc(assoc(k, 3));
            
            V = V + get_square(xi, yi, xj, yj, xk, yk)
        end
        
        fe(2*a+1) = fe(2*a+1) - d1*9.8*V*1/3;
     end
    for i=1:size(ige2,2)
        tmp = [];
        a = ige2(i);
        for j=1:els
            it = assoc(j, 1);
            jt = assoc(j, 2);
            kt = assoc(j, 3);
            if(a == it || a==jt || a==kt)
                tmp = [tmp j];
            end
        end
        V = 0;
        for k=1:size(tmp,2)
            xi = xc(assoc(k, 1));
            yi = yc(assoc(k, 1));

            xj = xc(assoc(k, 2));
            yj = yc(assoc(k, 2));

            xk = xc(assoc(k, 3));
            yk = yc(assoc(k, 3));
            
            V = V + get_square(xi, yi, xj, yj, xk, yk)
        end
        
        fe(2*a+1) = fe(2*a+1) - d1*9.8*V*1/3;
    end
    fe1=fe;
end

