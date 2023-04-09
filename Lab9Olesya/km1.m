function km1 = km1(ige1m, assoc, nds, xnds, ynds, E1,v1)
    kol_els_1 = size(ige1m,2)
    K = zeros(2*nds, 2*nds);
    for i=1:kol_els_1
        xi = xnds(assoc(ige1m(i),1));
        yi = ynds(assoc(ige1m(i),1));

        xj = xnds(assoc(ige1m(i),2));
        yj = ynds(assoc(ige1m(i),2));

        xk = xnds(assoc(ige1m(i),3));
        yk = ynds(assoc(ige1m(i),3));
    
        D = E1 * (1 - v1) / ((1 + v1) * (1 - 2 * v1)) * [[1, v1 / (1 - v1), 0]; [v1 / (1 - v1), 1, 0];[0, 0, (1 - 2 * v1) / (2 * (1 - v1))]];
       
        der = [[-1,0,1];[-1,1,0]];
        J = der*[[xi, yi]; [xj, yj];[xk, yk]];
        
        B=zeros(3,6);
        invJ = inv(J);
        
        tmp = invJ*der(:, 1);
        B(1, 1) = tmp(1);
        B(3, 1) = tmp(2);
        B(2, 2) = tmp(2);
        B(3, 2) = tmp(1);
        
        tmp = invJ * der(:, 2);
        
        B(1, 3) = tmp(1);
        B(3, 3) = tmp(2);
        B(2, 4) = tmp(2);
        B(3, 4) = tmp(1);
        
        tmp = invJ * der(:,3);
        
        B(1, 5) = tmp(1);
        B(3, 5) = tmp(2);
        B(2, 6) = tmp(2);
        B(3, 6) = tmp(1);
        
        S = get_square(xi, yi, xj, yj, xk, yk);

        ke = B'*D*B*S     
        
        A = zeros(6, 2*nds);

        it = assoc(ige1m(i), 1);
        jt = assoc(ige1m(i), 2);
        kt = assoc(ige1m(i), 3);

        A(1, 2*it-1) = 1;
        A(2, 2*it) = 1;
        A(3, 2 * jt-1) = 1;
        A(4, 2 * jt) = 1;
        A(5, 2 * kt-1) = 1;
        A(6, 2 * kt ) = 1;
        K = K + A'*ke*A;
    end
    km1 = K;
end

