function  invphi = mainpass(matrices,add,fdd,gdd,D2)
%% Initial store for body inverse inertia                          

        invphi.PHI_STORE11{1,1} = matrices.body1_zeta11;
        invphi.PHI_STORE22{1,1} = matrices.body1_zeta22;
        invphi.PHI_STORE12{1,1} = matrices.body1_zeta12;
        invphi.PHI_STORE21{1,1} = matrices.body1_zeta21;
        invphi.PHI_STORE13{1,1} = matrices.body1_zeta13;
        invphi.PHI_STORE23{1,1} = matrices.body1_zeta23;
        
        invphi.PHI_STORE11{1,2} = matrices.body2_zeta11;
        invphi.PHI_STORE22{1,2} = matrices.body2_zeta22;
        invphi.PHI_STORE12{1,2} = matrices.body2_zeta12;
        invphi.PHI_STORE21{1,2} = matrices.body2_zeta21;
        invphi.PHI_STORE13{1,2} = matrices.body2_zeta13;
        invphi.PHI_STORE23{1,2} = matrices.body2_zeta23;


        
        invphi.Xhat{2,1} = D2*inv(D2'*(invphi.PHI_STORE22{1,1}+invphi.PHI_STORE11{1,2})*D2)*D2';
        invphi.Y{2,1}    = (invphi.PHI_STORE23{1,1}-invphi.PHI_STORE13{1,2}+[add(2)+fdd;0;0]);
        
        invphi.PHI_STORE11{2,1} = invphi.PHI_STORE11{1,1}-invphi.PHI_STORE12{1,1}*invphi.Xhat{2,1}*invphi.PHI_STORE21{1,1};
        invphi.PHI_STORE22{2,1} = invphi.PHI_STORE22{1,2}-invphi.PHI_STORE21{1,2}*invphi.Xhat{2,1}*invphi.PHI_STORE12{1,2};
        
        invphi.PHI_STORE12{2,1} = invphi.PHI_STORE12{1,1}*invphi.Xhat{2,1}*invphi.PHI_STORE12{1,2};
        invphi.PHI_STORE21{2,1} = invphi.PHI_STORE21{1,2}*invphi.Xhat{2,1}*invphi.PHI_STORE21{1,1};
        
        invphi.PHI_STORE13{2,1} = invphi.PHI_STORE13{1,1}-invphi.PHI_STORE12{1,1}*invphi.Xhat{2,1}*invphi.Y{2,1};
        invphi.PHI_STORE23{2,1} = invphi.PHI_STORE23{1,2}+invphi.PHI_STORE21{1,2}*invphi.Xhat{2,1}*invphi.Y{2,1};
