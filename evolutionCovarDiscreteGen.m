function [G11, G22, G12, betaVector, GAngle, GEpsilon,numMutOverTime,regressionOverTime]=evolutionCovarDiscreteGen(w11, w22, sCor, tT, aT, steps,N)

% After Jones 2003
% STABILITY OF THE G-MATRIX IN A POPULATION EXPERIENCE PLEITROPIC MUTAITON,
% STABILIZING SELECTION, AND GENETIC DRIFT
isScript=false;
verbose=false;
plotEllipse=false;
environment=true;
randomInit=false;
if(isScript)
    tT=0; % transient time
    aT=1000; % time for analysis
    steps=10;
    
    N=256; % num individuals
    sCor=0; % r_w -- selection correlation
    
end
mCor=0; % r_mu -- mutation correlation
n=50; % number of loci


mutProb=0.0002; % probability of mutation per loci



mVar1=0.05; % alpha_1^2 -- the mutation variance for trait 1
mVar2=0.05; % alpha_2^2

mCoVar=mCor*sqrt(mVar1)*sqrt(mVar2);
mutationMatrix=[mVar1, mCoVar; mCoVar, mVar2]; % the mutation matrix for the bivariate normal
mutMean=zeros(1,2); % the mutation mean -- mutation bias can be implemented here

if(isScript)
    w11=49; % selection variance for trait 1
    w22=9; % selection variance for trait 2
end
w12=sCor*sqrt(w11)*sqrt(w22); % selection covariance
selectionMatrix=[w11, w12; w12, w22]; % the matrix omega (w) in the paper
% describes the curvature and orientation of the selection surface
% progeny=zeros(n,2,2);

numMutOverTime=zeros(1,aT); % number of mutants over time
regressionOverTime=zeros(2,aT); % regression of parent to child trait over time
%########################################################################%
popMatrix=zeros(N,n,2, 2);
%     Tracks the effects of the ith loci on the jth individual of the kth
%     chromosomal set of the lth trait in popMatrix(j, i, k, l)
if(randomInit)
    for i=1:N
        for j=1:n
            popMatrix(i, j, 1, :)=mvnrnd(mutMean, mutationMatrix);
            popMatrix(i, j, 2, :)=mvnrnd(mutMean, mutationMatrix);
        end
    end
    
end

popTrait=zeros(N, 2); % population trait values
if(randomInit)
    for i=1:N
        popTrait(i,1)=sum(popMatrix(i, :, 1));%+randn; % randn is the environment
        popTrait(i,2)=sum(popMatrix(i, :, 2));%+randn;
    end
end
%########################################################################%


fitness=zeros(1,N); % track every agent's fitness


optTraitVal=zeros(1,2); % theta, the optimal trait value for the 2 traits



for agent=1:N
    % calculate fitness
    
    fitness(agent)=calcFitness(popTrait(agent,:), ...
        optTraitVal, selectionMatrix);
end
% avgFitness=mean(fitness);


GAngle=zeros(1,aT); % the angle of the eigenvector of G from the x-axis
GEpsilon=zeros(1,aT); % the shape parameter, 1/eccentricity
G11=zeros(1,aT); % variance in trait 1, G11
G22=zeros(1,aT); % variance in trait 2, G22
G12=zeros(1,aT); % covariance in trait 1 and 2, G12
% G21=zeros(1,aT);

betaVector=zeros(2,aT);
% measures the covariance between the trait values and fitness over time

parentChildCor=zeros(5*N, 2, 2);
% this keep track of parent-child correlation for calculating the G matrix
% tracks the ith birth event's jth trait in (i,j,1) for parent, (i,j,2)
% for child
% numChild=1;


if(plotEllipse)
    ellipseFigure=figure;
end

trackParentChild=false;

progenyPopMatrix=zeros(2*N,n,2, 2); % the matrix of all progeny
progenyTraitMatrix=zeros(2*N, 2); % traits of all progeny
progenyFitness=zeros(1,2*N);

for t=1:tT+aT
    if(t>tT) % last time step
        trackParentChild=true;
    end
    numMut=0;
    
    numChild=0;
    % tracks the number of the child born, index for the parentChildCor
    
    mateOrder=randperm(N); % if all are to mate
    
    
    
    for i=1:N/2
        
        
        
        
        firstAgent=mateOrder(2*i-1);
        secondAgent=mateOrder(2*i); % the two monogamous mates
        
        %         firstAgent=ceil(rand*N);
        %         secondAgent=ceil(rand*N);
        %         while (secondAgent==firstAgent)
        %             secondAgent=ceil(rand*N);
        %         end
        
        %         fit1=fitness(firstAgent);
        %         fit2=fitness(secondAgent);
        if(trackParentChild)
            firstParentTrait1=popTrait(firstAgent, 1);
            secondParentTrait1=popTrait(secondAgent, 1);
            firstParentTrait2=popTrait(firstAgent, 2);
            secondParentTrait2=popTrait(secondAgent, 2);
        end
        %         meanFit=(fit1+fit2)/2;
        %         relFit=meanFit;
        %         relFit=meanFit/mean(fitness);
        
        
        if(verbose)
            disp(['first agent: ' num2str(firstAgent)]);
            disp(['first fitness: ' num2str(fit1)]);
            disp(['second agent: ' num2str(secondAgent)]);
            disp(['second fitness: ' num2str(fit2)]);
            disp(['mean fitness: ' num2str(mean(fitness))]);
            disp(['relative fitness: ' num2str(relFit)]);
            %             disp(['number of progeny: ' num2str(numProgeny)]);
        end
        
        for ii=1:4 % they have 4 children, as per Jones
            %             if(rand<relFit) % then the child is born, else not viable
            
            %             toMut=rand(2,n)<mutProb; % mutation on child chromosome set
            
            
            
            chromChoice1=(rand(1,n)<0.5)+1; % which allele from parent 1 -- half from chromosome 1, half from chromosome 2
            chromChoice2=(rand(1,n)<0.5)+1; % which allele from parent 2
            
            allelesParent1=sub2ind(size(popMatrix), ones(1,2*n)*firstAgent, [1:n,1:n], [chromChoice1,chromChoice1], [ones(1,n), ones(1,n)*2]);
            allelesParent2=sub2ind(size(popMatrix), ones(1,2*n)*secondAgent, [1:n,1:n], [chromChoice2,chromChoice2], [ones(1,n), ones(1,n)*2]);
            
            progeny=zeros(n,2,2); % (allele, chromosome set, trait)
            
            progeny(:,1,:)=reshape(popMatrix(allelesParent1), n, 1, 2);
            progeny(:,2,:)=reshape(popMatrix(allelesParent2), n, 1, 2);
            
            
            
            
            %
            for iii=1:n
                %                 progeny(iii,1,:)=popMatrix(firstAgent, iii, chromChoice1,:);
                %                 progeny(iii,2,:)=popMatrix(secondAgent,iii, chromChoice2,:);
                %                 if(rand<0.5) % recombine at this loci?
                %                     temp=progeny(iii,1,:);
                %                     progeny(iii,1,:)=progeny(iii,2,:); % swap loci of the two chromosomes
                %                     progeny(iii,2,:)=temp;
                %                 end
                %
                %                     progeny(iii,1,:)=reshape(popMatrix(firstAgent,iii,(rand<0.5)+1,:), 1, 1, 2); % allele on first chromosome from parent 1, but random parent 1 chromosome
                %                     progeny(iii,2,:)=reshape(popMatrix(secondAgent,iii,(rand<0.5)+1, :), 1, 1, 2);
                for iiii=1:2
                    if(rand<mutProb)
                        numMut=numMut+1;
                        mutAmount=mvnrnd(mutMean, mutationMatrix);
                        progeny(iii,iiii,:)=progeny(iii,iiii,:)+reshape(mutAmount, 1, 1, 2);
                    end
                end
                %
            end
            
            
            
            progenyTrait=zeros(1,2);
            progenyTrait(1)=sum(progeny(:,1, 1)+progeny(:, 2, 1));
            progenyTrait(2)=sum(progeny(:,1, 2)+progeny(:, 2, 2));
            if(environment)
                progenyTrait(1)=progenyTrait(1)+randn;
                progenyTrait(2)=progenyTrait(2)+randn;
            end
            
            
            
            
            if(verbose)
                disp(['alleles to mutate: ' num2str(toMut)]);
                disp(['Choice of alleles: ' num2str(alleleChoice)]);
                disp(['first parent allele effects on trait 1: ' num2str(popMatrix(firstAgent,:,1))]);
                disp(['first parent allele effects on trait 2: ' num2str(popMatrix(firstAgent,:,2))]);
                disp(['second parent allele effects on trait 1: ' num2str(popMatrix(secondAgent,:,1))]);
                disp(['second parent allele effects on trait 1: ' num2str(popMatrix(secondAgent,:,2))]);
                disp(['first parent traits: ' num2str([firstParentTrait1, firstParentTrait2])]);
                disp(['second parent traits: ' num2str([secondParentTrait1, secondParentTrait2])]);
                disp(['Child allele effects on trait 1: ' num2str(progeny(:,1)')]);
                disp(['Child allele effects on trait 2: ' num2str(progeny(:,2)')]);
                disp(['Child allele trait 1: ' num2str(sum(progeny(:,1)))]);
                disp(['Child allele trait 2: ' num2str(sum(progeny(:,2)))]);
                %                 disp(['Parent-child correlation matrix: ' num2str(parentChildCor(numChild, :))]);
                
                
            end
            thisProgenyFitness=calcFitness(progenyTrait, optTraitVal, selectionMatrix);
            if(rand<thisProgenyFitness)
                
                numChild=numChild+1;
                progenyFitness(numChild)=thisProgenyFitness;
                progenyPopMatrix(numChild, :,:,:)=progeny;
                progenyTraitMatrix(numChild,:)=progenyTrait;
            end
            if(trackParentChild && numChild>0)
                parentChildCor(numChild, 1, 1)=(firstParentTrait1+secondParentTrait1)/2; % midpoint value = parent trait 1
                parentChildCor(numChild, 2, 1)=(firstParentTrait2+secondParentTrait2)/2; % midpoint value = parent trait 2
                
                parentChildCor(numChild, 1,2)=progenyTrait(1); % child trait 1
                parentChildCor(numChild, 2,2)=progenyTrait(2); % child trait 2
                
            end
        end
        %         end
        
        
        
        
        
    end
    if(mod(t,steps)==0)
        disp(['time: ' num2str(t), ' number of progeny: ' num2str(numChild), ' number of mutants: ' num2str(numMut)]);
    end
    
    temp2=randperm(numChild);
    survivors=temp2(1:N); % this will quit with error if numChild<survivors
    popMatrix=progenyPopMatrix(survivors,:,:,:);
    popTrait=progenyTraitMatrix(survivors,:);
%     for numSurvivor=1:N
    fitness=progenyFitness(survivors);
%     end
    
    if(t>tT)
        %
        %         g=cov([parentChildCor(1:numChild, 1, 1), parentChildCor(1:numChild, 1, 2)]);
        %         G11(t-tT)=g(1, 2);
        %         g=cov([parentChildCor(1:numChild, 2, 1), parentChildCor(1:numChild, 2, 2)]);
        %         G22(t-tT)=g(1,2);
        %         g=cov([parentChildCor(1:numChild, 1, 1), parentChildCor(1:numChild, 2, 2)]);
        %         G12(t-tT)=g(1,2);
        %         g=cov([parentChildCor(1:numChild, 2, 1), parentChildCor(1:numChild, 1, 2)]);
        %         G21(t-tT)=g(1,2);
        %         G=[G11(t-tT), G12(t-tT); G21(t-tT), G22(t-tT)];
        %
        numMutOverTime(t-tT)=numMut;
        p1= polyfit(parentChildCor(1:numChild-2, 1, 1), parentChildCor(1:numChild-2, 1, 2),1);
        p2= polyfit(parentChildCor(1:numChild-2, 2, 1), parentChildCor(1:numChild-2, 2, 2),1);
        
        regressionOverTime(:,t-tT)=[p1(1); p2(1)];
        beta1=cov(fitness, popTrait(:,1));
        beta2=cov(fitness, popTrait(:,2));
        betaVector(:,t-tT)=[beta1(1, 2), beta2(1, 2)]';
        
        breedingValues=reshape(sum(sum(popMatrix, 3), 2), N, 2);
        % trait value is the sum across the two chromosomal sets, then across each allele
        
        covarianceMatrix=cov(breedingValues);
        G11(t-tT)=covarianceMatrix(1,1);
        G22(t-tT)=covarianceMatrix(2,2);
        G12(t-tT)=covarianceMatrix(1,2);
        G=[G11(t-tT), G12(t-tT); G12(t-tT), G22(t-tT)];
        
        
        [eigVectors, eigValues]=eig(G);
        eigVal1=abs(eigValues(1,1));
        eigVal2=abs(eigValues(2,2));
        [bigEigVal, bigIndex]=max([eigVal1, eigVal2]);
        [smallEigVal, smallIndex]=min([eigVal1, eigVal2]);
        GEpsilon(t-tT)=smallEigVal/bigEigVal;
        bigVector=eigVectors(:, bigIndex);
        smallVector=eigVectors(:,smallIndex);
        GAngle(t-tT)=radtodeg(atan2(real(bigVector(2)), real(bigVector(1))));
        if(GAngle(t-tT)<-90)
            GAngle(t-tT)=GAngle(t-tT)+180;
        elseif(GAngle(t-tT)>90)
            GAngle(t-tT)=GAngle(t-tT)-180;
        end
        if(mod(t,steps)==0)
            disp([num2str(t) ' -- GAngle: ' num2str(GAngle(t-tT)), ' GEpsilon: ' num2str(GEpsilon(t-tT)), ' GSize: ' num2str(G11(t-tT)+G22(t-tT)), ' numMut ' num2str(numMut)]);
        end
        if(plotEllipse)
            p=calculateEllipse(0,0, smallEigVal, bigEigVal, GAngle(t-tT));
            figure(ellipseFigure);
            
            plot(p(:,1), p(:,2), '.-');
            axis([-2 2 -2 2]);
        end
        
    end
    
end

