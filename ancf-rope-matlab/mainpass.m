function  invphi = mainpass(matrices,ne)

%% Initial store for body inverse inertia                          
        for i = 1:ne
        invphi.PHI_STORE11{1,i} = matrices.body_zeta11{i};
        invphi.PHI_STORE22{1,i} = matrices.body_zeta22{i};
        invphi.PHI_STORE12{1,i} = matrices.body_zeta12{i};
        invphi.PHI_STORE21{1,i} = matrices.body_zeta21{i};
        invphi.PHI_STORE13{1,i} = matrices.body_zeta13{i};
        invphi.PHI_STORE23{1,i} = matrices.body_zeta23{i};
        end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        for pass = 2:ne    
        invphi.Xhat{pass,1} = (invphi.PHI_STORE22{pass-1,1}+invphi.PHI_STORE11{1,pass});
        invphi.W{pass,1}    = inv(invphi.Xhat{pass,1});
        invphi.Y{pass,1}    = invphi.W{pass,1}*(invphi.PHI_STORE23{pass-1,1}-invphi.PHI_STORE13{1,pass});
        
        invphi.PHI_STORE11{pass,1} = invphi.PHI_STORE11{pass-1,1}-invphi.PHI_STORE12{pass-1,1}*invphi.W{pass,1}*invphi.PHI_STORE21{pass-1,1};
        invphi.PHI_STORE22{pass,1} = invphi.PHI_STORE22{1,pass}-invphi.PHI_STORE21{1,pass}*invphi.W{pass,1}*invphi.PHI_STORE12{1,pass};
        
        invphi.PHI_STORE12{pass,1} = invphi.PHI_STORE12{pass-1,1}*invphi.W{pass,1}*invphi.PHI_STORE12{1,pass};
        invphi.PHI_STORE21{pass,1} = invphi.PHI_STORE21{1,pass}*invphi.W{pass,1}*invphi.PHI_STORE21{pass-1,1};
        
        invphi.PHI_STORE13{pass,1} = invphi.PHI_STORE13{pass-1,1}-invphi.PHI_STORE12{pass-1,1}*invphi.Y{pass,1};
        invphi.PHI_STORE23{pass,1} = invphi.PHI_STORE23{1,pass}+invphi.PHI_STORE21{1,pass}*invphi.Y{pass,1};
        end
