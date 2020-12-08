function [Trajectory,timeVector] = PlasmidSpreadFunction2d_PaintTree(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,simID,ChangeTime,InitialRatios)

%SharingNess=[0,0.2,0.001];
q=1;

timeVector=zeros(1,8);
TSTART=tic();

N=2^Q;
Jumps=[1,-1,N,-N];

ChromosoneType= rand(N^2,1);

ChromosoneType= 1*(ChromosoneType<InitialRatios(1))+2*1*(ChromosoneType>1-InitialRatios(2));

PlasmidType=randi(k,N^2,q);
PlasmidMultiplier= k.^((0:q))';

ColorWheel= [eye(3);ones(3)-eye(3)];
ChromosoneMultiplier=[0,1,0.5];

SavePaints={};
SavePaintsCount=1;

paintTimes=[0,500,1000,2000,5000,inf];
paintIndex=1;

Filter=zeros(k^q,k*q);
for(qqq=1:q)
    Filter(:,(1:k)+k*(qqq-1))= kron(kron(ones(k^(qqq-1),1) ,eye(k,k) ), ones(k^(q-qqq),1));
end
Filter=Filter==1;


Environment= mod(mod( floor((1:N^2)/N),2)+(1:N^2),2)' ;%%Deafault to checkerboard
if(EnviroShape==0)
    Environment= 0*Environment; 
elseif(EnviroShape==1)
    Environment= 1*((1:N^2)>N/2)'; %A basic 2 environs set up. Other options possible.
elseif(EnviroShape==3)
    Environment= (rand(1,N^2)>0.5)'; %Random Environment. Yay!
end
    
    
fitnessBox= [baseFit+FitDeltaCombo*rand(k^q,1)+Filter*FitDeltaLin*rand(k*q,1);baseFit+FitDeltaCombo*rand(k^q,1)+Filter*FitDeltaLin*rand(k*q,1)];
DeathRate = [baseDeath+FitDeltaDeath*rand(k^q,1);baseDeath+FitDeltaDeath*rand(k^q,1)];
%DeathRate=[DeathRate;DeathRate];



t=0;


Tstop=20000;
States=zeros(length(fitnessBox)+7,Tstop+1);
bookmark=0;

toRefresh=(1:N^2);

FullIndex=([PlasmidType(toRefresh,:)-1,Environment(toRefresh)]*PlasmidMultiplier+1).*(ChromosoneType(toRefresh)>0)+1*(ChromosoneType(toRefresh)==0);
Rates=[fitnessBox(FullIndex),DeathRate(FullIndex),mutationRate*ones(length(FullIndex),1),SharingNess(ChromosoneType(toRefresh)+1)'*ones(1,q)]';
Rates=Rates.*(ChromosoneType(:)>0)';

Weights= sum(Rates);

treeWeights=zeros(1,length(Weights));

treeWeights(end)=sum(Weights);

for( qqq=(1:2*Q) )
   for(aaa=1:(2^(qqq-1)) )
       treeWeights(2^(qqq-1)-1 + aaa)=sum(Weights( (2*(2^(2*Q-qqq))*(aaa-1) +1):(2*(2^(2*Q-qqq))*(aaa-0.5)) ));
   end
end


CumSumRates=cumsum(Rates(:));

toRefresh=[];

numEvents=0;
tic()
expRoller=exprnd(1,10003,1);
PickRoller=rand(10003,1);
ActionRoller=rand(10003,1);
timeVector(2)=timeVector(2)+toc();
RollerDex=1;

deathCount=0;
birthCount=0;
transferCount=0;
mutantCount=0;

fitnesses= zeros(k,Tstop+1);
ChromDensity= zeros(2,Tstop+1);
PlasmidDensity=zeros(k,Tstop+1);

SavePaints={}

        R=(ColorWheel(PlasmidType+1,1).*ChromosoneMultiplier(ChromosoneType+1)')>0;
        G=(ColorWheel(PlasmidType+1,2).*ChromosoneMultiplier(ChromosoneType+1)')>0;
        B=(ColorWheel(PlasmidType+1,3).*ChromosoneMultiplier(ChromosoneType+1)')>0;
        A= ChromosoneType==2;
        
        image=false(N,N,4); %initialize
        image(:,:,1)=reshape(R,N,N);
        image(:,:,2)=reshape(G,N,N);
        image(:,:,3)=reshape(B,N,N);
        image(:,:,4)=reshape(A,N,N);

        SavePaints{SavePaintsCount}=image;
        SavePaintsCount=1+SavePaintsCount;        

        
while(t<Tstop)

    NextTime= expRoller(RollerDex)./treeWeights(end);
    
    target= PickRoller(RollerDex)*treeWeights(end); 
    
        Index=1;
        while(Index<length(treeWeights))
            if(target>treeWeights(Index))
                target=target-treeWeights(Index);
                Index=Index*2+1;
            else
                Index=Index*2;
            end    
        end

        NextIndex=Index-length(treeWeights)+1;
        targ=ActionRoller(RollerDex)*sum(Rates(:,NextIndex));
        ActionPick=sum( targ>cumsum(Rates(:,NextIndex)));
     RollerDex=RollerDex+1;
    
    if(t+NextTime>paintTimes(paintIndex))
       
        R=ColorWheel(PlasmidType+1,1).*ChromosoneMultiplier(ChromosoneType+1)';
        G=ColorWheel(PlasmidType+1,2).*ChromosoneMultiplier(ChromosoneType+1)';
        B=ColorWheel(PlasmidType+1,3).*ChromosoneMultiplier(ChromosoneType+1)';
        
         image=zeros(N,N,3); %initialize
        image(:,:,1)=reshape(R,N,N);   %Red (dark red)
        image(:,:,2)=reshape(G,N,N);
        image(:,:,3)=reshape(B,N,N);
        figure(1)
        subplot(1,length(paintTimes)-1,paintIndex);
        imshow(image);
        
        paintIndex=paintIndex+1;
       
    end
     
 if(t+NextTime>bookmark)
     
    t=bookmark;
    NextTime=0;
    bookmark=bookmark+1;
    
    ActionPick=-ActionPick-0.1;
    
       toRefresh=1:N^2;
       
    FullIndex=([PlasmidType(toRefresh,:)-1,Environment(toRefresh)]*PlasmidMultiplier+1).*(ChromosoneType(toRefresh)>0)+1*(ChromosoneType(toRefresh)==0);
    
    for(qqq=1:length(fitnessBox))
        States(qqq,bookmark)= sum(FullIndex==qqq&(ChromosoneType>=1))/N^2;
    end
    
    States(length(fitnessBox)+1,bookmark) = mean(fitnessBox(FullIndex(ChromosoneType>0) ));
    States(length(fitnessBox)+2,bookmark) = mean(DeathRate(FullIndex(ChromosoneType>0) ));
    States(length(fitnessBox)+3,bookmark) = mean(fitnessBox(FullIndex(ChromosoneType>0) )./DeathRate(FullIndex(ChromosoneType>0) ));
    States(length(fitnessBox)+4,bookmark) = mean(ChromosoneType==1);
    States(length(fitnessBox)+5,bookmark) = mean(ChromosoneType==2);

        ChromDensity(1,bookmark)= mean(ChromosoneType==1);
        ChromDensity(2,bookmark)= mean(ChromosoneType==2);
        
        for(kkk=1:k)
            PlasmidDensity(kkk,bookmark)= mean(PlasmidType==kkk);
            fitnesses(kkk,bookmark)=fitnessBox(kkk)/DeathRate(kkk);
        end
%    StatesSorted(:,bookmark)= [FilterA;FilterB;FilterC]*States(1:64,bookmark);

    if(States(length(fitnessBox)+4,bookmark)==0 || States(length(fitnessBox)+5,bookmark)==0)
            States = States(:,1:bookmark);
 %           StatesSorted=StatesSorted(:,1:bookmark);
            break;
    end
   timeVector(3)=timeVector(3)+toc();
    
    if(mod(bookmark,ChangeTime)<1)
        
        R=(ColorWheel(PlasmidType+1,1).*ChromosoneMultiplier(ChromosoneType+1)')>0;
        G=(ColorWheel(PlasmidType+1,2).*ChromosoneMultiplier(ChromosoneType+1)')>0;
        B=(ColorWheel(PlasmidType+1,3).*ChromosoneMultiplier(ChromosoneType+1)')>0;
        A= ChromosoneType==2;
        
        image=false(N,N,4); %initialize
        image(:,:,1)=reshape(R,N,N);
        image(:,:,2)=reshape(G,N,N);
        image(:,:,3)=reshape(B,N,N);
        image(:,:,4)=reshape(A,N,N);

        SavePaints{SavePaintsCount}=image;
        SavePaintsCount=1+SavePaintsCount;        
        
            [numEvents,bookmark,simID]
       if(rand()>0.5)
            fitnessBox(1:PlasmidMultiplier(end))=  [baseFit+FitDeltaCombo*rand(k^q,1)+Filter*FitDeltaLin*rand(k*q,1)];
            DeathRate(1:PlasmidMultiplier(end)) = [baseDeath+FitDeltaDeath*rand(k^q,1)];    
       else
            fitnessBox((PlasmidMultiplier(end)+1):end)= [baseFit+FitDeltaCombo*rand(k^q,1)+Filter*FitDeltaLin*rand(k*q,1)]';
            DeathRate((PlasmidMultiplier(end)+1):end) = [baseDeath+FitDeltaDeath*rand(k^q,1)];
       end
       
        toRefresh=1:N^2;
        FullIndex=([PlasmidType(toRefresh,:)-1,Environment(toRefresh,1)]*PlasmidMultiplier+1).*(ChromosoneType(toRefresh)>0)+1*(ChromosoneType(toRefresh)==0);
        Rates=[fitnessBox(FullIndex),DeathRate(FullIndex),mutationRate*ones(length(FullIndex),1),SharingNess(ChromosoneType(toRefresh)+1)'*ones(1,q)]';
        Rates=Rates.*(ChromosoneType(:)>0)';
        Weights= sum(Rates);
        treeWeights=zeros(1,length(Weights));
        treeWeights(end)=sum(Weights);

        for( qqq=(1:2*Q) )
            for(aaa=1:(2^(qqq-1)) )
                treeWeights(2^(qqq-1)-1 + aaa)=sum(Weights( (2*(2^(2*Q-qqq))*(aaa-1) +1):(2*(2^(2*Q-qqq))*(aaa-0.5)) ));
            end
        end
        
    end
    
    toRefresh=ones(0,1);
    
elseif(ActionPick==1) %%Death
    ChromosoneType(NextIndex)=0;
    PlasmidType(NextIndex,:)=0;
    toRefresh=NextIndex;
    deathCount=deathCount+1;
elseif(ActionPick==0) %%Birth
    Reciever= mod(NextIndex+Jumps(randi(length(Jumps)))-1,N^2)+1;
    if(ChromosoneType(Reciever)==0)
        ChromosoneType(Reciever)=ChromosoneType(NextIndex);
        PlasmidType(Reciever,:)=PlasmidType(NextIndex,:);
        toRefresh=[Reciever];
        birthCount=birthCount+1;
    else
        toRefresh=ones(0,1);
    end
elseif(ActionPick==2) %%Mutation
   Reciever=randi(q);
   PlasmidType(NextIndex,Reciever)=randi(k);
   toRefresh=NextIndex; 
    mutantCount=mutantCount+1;
else
    Reciever= mod(NextIndex+Jumps(randi(length(Jumps)))-1,N^2)+1;
    if(ChromosoneType(Reciever)>0)
        PlasmidType(Reciever,ActionPick-2)=PlasmidType(NextIndex,ActionPick-2);
        toRefresh=Reciever;
        transferCount=transferCount+1;
    else
        toRefresh=ones(0,1);
    end
 end

t=t+NextTime;
numEvents=numEvents+1;

if(~isempty(toRefresh))

    oldRate=sum(Rates(:,toRefresh));
   if(ChromosoneType(toRefresh)==0)
        Rates(:,toRefresh)=0;       
   else
       FullIndex=([PlasmidType(toRefresh,:)-1,Environment(toRefresh,1)]*PlasmidMultiplier+1).*(ChromosoneType(toRefresh)>0)+1*(ChromosoneType(toRefresh)==0);
        Rates(:,toRefresh)=[fitnessBox(FullIndex),DeathRate(FullIndex),mutationRate*ones(length(FullIndex),1),SharingNess(ChromosoneType(toRefresh)+1)'*ones(1,q)]';       
   end
    RateDiff=sum(Rates(:,toRefresh))-oldRate;

    treeWeights(end)=treeWeights(end)+RateDiff;
    DivisorThing=toRefresh+length(treeWeights)-1;
    while(DivisorThing>1)
        DivisorThing=DivisorThing/2;
        if(DivisorThing==floor(DivisorThing))
            treeWeights(DivisorThing)=treeWeights(DivisorThing)+RateDiff;
        end
        DivisorThing=floor(DivisorThing);
    end
    
end   


%      %%%%%%%%%%%%% TEST BINARY TREE ACCURACY
%      
%         FullIndex=([PlasmidType(:,:)-1,Environment(:,1)]*PlasmidMultiplier+1).*(ChromosoneType(:)>0)+1*(ChromosoneType(:)==0);
%         TestRates=[fitnessBox(FullIndex),DeathRate(FullIndex),mutationRate*ones(length(FullIndex),1),SharingNess(ChromosoneType+1)'*ones(1,q)]';
%         TestRates=TestRates.*(ChromosoneType(:)>0)';
%         Weights= sum(TestRates);
%         testTreeWeights=zeros(1,length(Weights));
%         testTreeWeights(end)=sum(Weights);
% 
%         for( qqq=(1:2*Q) )
%             for(aaa=1:(2^(qqq-1)) )
%                 testTreeWeights(2^(qqq-1)-1 + aaa)=sum(Weights( (2*(2^(2*Q-qqq))*(aaa-1) +1):(2*(2^(2*Q-qqq))*(aaa-0.5)) ));
%             end
%         end
% 
%         if( any(abs(testTreeWeights-treeWeights)>10^-8))
%             error('Tree is out of synch!');
%         end
    %%% END TEST BINARY TREE ACCURACY
        

if(RollerDex>length(PickRoller)-10)
    expRoller=exprnd(1,10002,1);
    PickRoller=rand(10002,1);
    RollerDex=1;
end

end
% 

figure(2)
subplot(3,1,1);
for(kkk=1:k)
    plot(fitnesses(kkk,1:bookmark)','LineWidth',2,'color',0.75*ColorWheel(kkk+1,:)); 
    hold on
end

ylabel('fitness \newline r(x)/d(x)');

subplot(3,1,2);
for(kkk=1:k)
    plot(PlasmidDensity(kkk,1:bookmark)','LineWidth',2,'color',0.75*ColorWheel(kkk+1,:)); 
    hold on
end

ylabel('Plasmids');

subplot(3,1,3);
plot(ChromDensity(1,1:bookmark)','LineWidth',2,'color',[1,0,0]);
hold on
plot(ChromDensity(2,1:bookmark)','LineWidth',2,'color',0.5*[1,0,0]);
ylabel('Chromosomes');
xlabel('time')

Trajectory=States(length(fitnessBox)+[4:5],:);
timeVector(1)= toc(TSTART);

end
% 
% target=15;
% 
% bip=SavePaints{target};
% image=zeros(64,64,3);
% image(:,:,1)= (1*bip(:,:,1)).*(1-0.4*bip(:,:,4));
% image(:,:,2)= (1*bip(:,:,2)).*(1-0.4*bip(:,:,4));
% image(:,:,3)= (1*bip(:,:,3)).*(1-0.4*bip(:,:,4));
% imshow(image)