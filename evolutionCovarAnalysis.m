finalG11=zeros(1,18);
finalG22=zeros(1,18);
finalG12=zeros(1,18);
finalGAngle=zeros(1,18);
finalGEpsilon=zeros(1,18);
delG11=zeros(1,18);
delG22=zeros(1,18);
delG12=zeros(1,18);
delGAngle=zeros(1,18);
delGEpsilon=zeros(1,18);


iter=0;
for k=1:3

    
    if(k==1)
        w11=9;
        w22=9;
    elseif(k==2)
        w11=49;
        w22=49;
    elseif(k==3)
        w11=49;
        w22=9;
    end
%     disp(['w11 ', num2str(w11) ' w22 ' num2str(w22)]);
    for sCor=[0, 0.25, 0.5, 0.75, 0.85, 0.9]

        iter=iter+1;
        load(['GM_DiscreteGen_', num2str(w11), '_', num2str(w22), '_', num2str(sCor), '.mat']);
        finalG11(iter)=mean(mean(G11_it));
        finalG22(iter)=mean(mean(G22_it));
        finalG12(iter)=mean(mean(G12_it));
        finalGAngle(iter)=mean(mean(GAngle_it));
        finalGEpsilon(iter)=mean(mean(GEpsilon_it));
        delG11(iter)=mean(mean(abs((G11_it(:,1:end-1)-G11_it(:,2:end))./G11_it(:,2:end))));
        delG22(iter)=mean(mean(abs((G22_it(:,1:end-1)-G22_it(:,2:end))./G22_it(:,2:end))));
        delG12(iter)=mean(mean(abs((G12_it(:,1:end-1)-G12_it(:,2:end)))));
        delGAngle(iter)=mean(mean(abs((GAngle_it(:,1:end-1)-GAngle_it(:,2:end)))));
        delGEpsilon(iter)=mean(mean(abs((GEpsilon_it(:,1:end-1)-GEpsilon_it(:,2:end))./GEpsilon_it(:,2:end))));
    end  
end
results=[finalG11',finalG22',finalG12',finalGAngle',finalGEpsilon',delG11',delG22',delG12',delGAngle',delGEpsilon'];
disp(results);
