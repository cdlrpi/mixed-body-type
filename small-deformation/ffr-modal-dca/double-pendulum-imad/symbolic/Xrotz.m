function  Rotz = Xrotz( theta )

c = cos(theta);
s = sin(theta);

Rotz = [  c  s  0  0  0  0 ;
         -s  c  0  0  0  0 ;
          0  0  1  0  0  0 ;
          0  0  0  c  s  0 ;
          0  0  0 -s  c  0 ;
          0  0  0  0  0  1];