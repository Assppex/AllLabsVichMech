function S = get_square(xi,yi,xj,yj,xk,yk)
    Sqm = zeros(3);
    Sqm(:, 3) = 1;
    Sqm(1, 1) = xi;
    Sqm(1, 2) = yi;
    Sqm(2, 1) = xj;
    Sqm(2, 2) = yj;
    Sqm(3, 1) = xk;
    Sqm(3, 2) = yk;
    S = 1 / 2 * abs(det(Sqm));
end

