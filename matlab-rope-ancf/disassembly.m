function  back = disassembly(invphi,ne)
% Disassembly Process
De     = [1, 0; 
          0, 1;
          0, 0;
          0, 0];
height = ne;
%% For the last assembly
back.f1{height,1}    = De*inv(De'*invphi.PHI_STORE11{height,1}*De)*(-De'*invphi.PHI_STORE13{height,1});
back.edd1{height,1}  = invphi.PHI_STORE11{height,1}*back.f1{height,1}+invphi.PHI_STORE13{height,1};
back.f2{height,1}    = zeros(4,1);
back.edd2{height,1}  = invphi.PHI_STORE21{height,1}*back.f1{height,1}+invphi.PHI_STORE23{height,1};


%% Disassembly Process
for pass = ne-1:-1:1

back.f1{pass,1}      =  back.f1{pass+1,1};
back.f2{1,pass+1}    =  back.f2{pass+1,1};
back.f1{1,pass+1}    =  invphi.W{pass+1,1}*(invphi.PHI_STORE21{pass,1}*back.f1{pass,1}-invphi.PHI_STORE12{1,pass+1}*back.f2{1,pass+1})+invphi.Y{pass+1,1};
back.f2{pass,1}      = -back.f1{1,pass+1};

back.edd1{pass,1}    =  back.edd1{pass+1,1};
back.edd2{1,pass+1}  =  back.edd2{pass+1,1};
back.edd2{pass,1}    =  invphi.PHI_STORE21{pass,1}*back.f1{pass,1}+invphi.PHI_STORE22{pass,1}*back.f2{pass,1}+invphi.PHI_STORE23{pass,1};
back.edd1{1,pass+1}  =  invphi.PHI_STORE11{1,pass+1}*back.f1{1,pass+1}+invphi.PHI_STORE12{1,pass+1}*back.f2{1,pass+1}+invphi.PHI_STORE13{1,pass+1};

end
