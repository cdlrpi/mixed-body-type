function matrices = mat(p, ne, iM, l, rho, A, Y, I)
                     
%% Body 1

for i = 1:ne
e  = p(:,i);
binv = bodyinertia(i, ne, e, iM, l, rho, A, Y, I);

matrices.body_zeta11{i} = binv.body_zeta11;
matrices.body_zeta12{i} = binv.body_zeta12;
matrices.body_zeta13{i} = binv.body_zeta13;
matrices.body_zeta21{i} = binv.body_zeta21;
matrices.body_zeta22{i} = binv.body_zeta22;
matrices.body_zeta23{i} = binv.body_zeta23;

end
