function  back = disassembly(invphi,add,fdd,gdd)
% Disassembly Process

back.A1{2,1}    = [add(1);0;0];
back.f1{2,1}    = inv(invphi.PHI_STORE11{2,1})*(back.A1{2,1}-invphi.PHI_STORE13{2,1});
back.f2{2,1}    = zeros(3,1);
back.A2{2,1}    = invphi.PHI_STORE21{2,1}*back.f1{2,1}+invphi.PHI_STORE22{2,1}*back.f2{2,1}+invphi.PHI_STORE23{2,1}+[gdd;0;0];

    
back.f1{1,1}      =  back.f1{2,1};
back.f2{1,2}      =  back.f2{2,1};
back.f1{1,2}      =  invphi.Xhat{2,1}*(invphi.PHI_STORE21{1,1}*back.f1{1,1}-invphi.PHI_STORE12{1,2}*back.f2{1,2}+invphi.Y{2,1});
back.f2{1,1}      = -back.f1{1,2};

back.A1{1,1}      =  back.A1{2,1};
back.A2{1,2}      =  back.A2{2,1};
back.A2{1,1}      =  invphi.PHI_STORE21{1,1}*back.f1{1,1}+invphi.PHI_STORE22{1,1}*back.f2{1,1}+invphi.PHI_STORE23{1,1}+[fdd;0;0];
back.A1{1,2}      =  invphi.PHI_STORE11{1,2}*back.f1{1,2}+invphi.PHI_STORE12{1,2}*back.f2{1,2}+invphi.PHI_STORE13{1,2};


