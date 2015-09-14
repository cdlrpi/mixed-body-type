function  pend = penmodelDCA(q1,q2)

pend.Xj0{1} = Xrotz(-q1);
pend.Xj0{2} = Xrotz(-(q1+q2));    

pend.X0j{1} = Xrotz(q1);
pend.X0j{2} = Xrotz(q1+q2);    