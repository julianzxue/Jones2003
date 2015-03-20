numRuns=20;
tT=10000; % transient time
aT=2000; % analysis time
steps=100;
N=256;
% [G11, G22, G12, betaVector, GAngle, GEpsilon]

G11_it=zeros(numRuns, aT);
G22_it=zeros(numRuns, aT);
G12_it=zeros(numRuns, aT);
betaVector_it=zeros(numRuns, 2, aT);
GAngle_it=zeros(numRuns, aT);
GEpsilon_it=zeros(numRuns, aT);

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
    disp(['w11 ', num2str(w11) ' w22 ' num2str(w22)]);
    for sCor=[0, 0.25, 0.5, 0.75, 0.85, 0.9]

        disp(['sCor: ', num2str(sCor)]);
        parfor it=1:numRuns
            [G11_it(it,:),G22_it(it,:),G12_it(it,:),betaVector_it(it,:,:),...
                GAngle_it(it,:),GEpsilon_it(it,:), numMutOverTime,...
                regressionOverTime ] = ...
                evolutionCovarDiscreteGen(w11, w22, sCor, tT, aT, steps,N);
            
        end
        disp(' ');
        disp(['Mean G11: ' num2str(mean(mean(G11_it)))]); 
        disp(['Mean G22: ' num2str(mean(mean(G22_it)))]);
        disp(['Mean G12: ' num2str(mean(mean(G12_it)))]);
        disp(['Mean GAngle: ' num2str(mean(mean(GAngle_it)))]);
        disp(['Mean GEpsilon: ' num2str(mean(mean(GEpsilon_it)))]);
        disp(' ');
        save(['GM_DiscreteGen_', num2str(w11), '_', num2str(w22), '_', num2str(sCor), '.mat']);
    end  
end
