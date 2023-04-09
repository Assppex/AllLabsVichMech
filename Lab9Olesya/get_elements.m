function [ige1mt, ige2mt, basemt, kol, kol1, kol2] = get_elements(xc, yc, assoc, els)
    ige1mt = [];
    ige2mt = [];
    basemt = [];
    kol = 0;
    kol1 = 0;
    kol2 = 0;
    for i=1: els
        xi = xc(assoc(i, 1));
        yi = yc(assoc(i, 1));

        xj = xc(assoc(i, 2));
        yj = yc(assoc(i, 2));

        xk = xc(assoc(i, 3));
        yk = yc(assoc(i, 3));

        xcm = 1/3*(xi + xj + xk);
        ycm = 1/3*(yi + yj + yk);

        if ((xcm < 9 / 5 * (ycm + 71 / 3) ) && ycm >= 0 && xcm >= 0)
            ige1mt = [ige1mt i];
            kol = kol + 1;
        end
        if ((xcm > 9 / 5 * (ycm + 71 / 3)) && ycm >= 0 && xcm <= 1029)
            ige2mt = [ige2mt i];
            kol1 = kol1 + 1;
        end
        if (ycm <= 0)
            basemt = [basemt i];
            kol2 = kol2 + 1;
        end
    end
end

