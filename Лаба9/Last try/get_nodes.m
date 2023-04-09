function [ige1,ige2,base,kol,kol1,kol2] = get_nodes(xnds,ynds,nds)
    kol = 0;
    kol1 = 0;
    kol2 = 0;
    ige1=[];
    ige2=[];
    base=[];
    for i=1:nds
        xi = xnds(i);
        yi = ynds(i);
        if ((xi < 9 / 5 * (yi + 71 / 3) || abs(yi - 5 / 9 * xi + 71 / 3) < 10^(-4)) && yi >= 0 && xi >= 0)
            ige1 = [ige1 i];
            kol = kol + 1;
        end
        if ((xi > 9 / 5 * (yi + 71 / 3) || abs(yi - 5 / 9 * xi + 71 / 3) < 10^(-4)) && yi >= 0 && xi <= 1029)
            ige2 = [ige2 i];
            kol1 = kol1 + 1;
        end
        if (yi <= 0)
            base = [base i];
            kol2 = kol2 + 1;
        end 
    end
end

