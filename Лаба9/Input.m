function [coordinates, N, nodes_1_1, nodes_1_2, nodes_2, BC_X_nodes, BC_Y_nodes, force_nodes, force_elemeTemp=readtable('xydata.txt');
xnds=Temp(:,2);
ynds=Temp(:,3);
xnds=table2array(xnds);
ynds=table2array(ynds);

Temp=readtable('assoc.txt'); 
a=Temp(:,2);
a=table2array(a);
b=Temp(:,3);
b=table2array(b);
c=Temp(:,4);
c=table2array(c);


nds=size(xnds,1);
els=size(Temp,1);nts] = Input()

end

