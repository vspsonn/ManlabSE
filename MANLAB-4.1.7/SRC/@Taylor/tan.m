function v = tan(u)
% Taylor/exp est la formule de reu.coefrrence pour l'exponetielle
global Ck

g=cos(u).*cos(u);
tan0=tan(u.value);

tu=u.coef*0; 
sto=1./g.value;
for k=1:Ck
    tu(:,:,k)=0;
    for j=1:k-1
        tu(:,:,k)=tu(:,:,k)+j*times(tu(:,:,j),g.coef(:,:,k-j));         
    end
    tu(:,:,k)=sto.*(u.coef(:,:,k)-tu(:,:,k)/k);
end
v = Taylor(u.order,tan0,tu);