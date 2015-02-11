function  edd = MainFunction(p, ne, iM, l, rho, A, Y, I)

matrices = mat(p, ne, iM, l, rho, A, Y, I);
invphi   = mainpass(matrices,ne);
back     = disassembly(invphi,ne);

edd = [];
for i = 1:ne
acc = [back.edd1{1,i};back.edd2{1,i}];
edd = [edd,acc];
end
