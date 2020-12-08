 

% % This is the central script used for itterating through agent based
% simulations. The prelude sets all parameter values, and then subsquent
% "windows" run 100 simulations for each of the various simulation
% conditions. 
% If you with to run only one experimental condition, then comment out all
% others.
% Runnking on a Standard "Thinkpad 13", each of these conditions takes 1-2
% hours to process. Running through all of them can take 2-3 days.
%
% If you want to know what is going on in each window then look for the
% "comment{?}= " command- the comment will tell you what experimental
% condition you are looking at.

%
%

close all;
Q=6;  %% N=64
k=4;
q=1;
mutationRate=10^-5;
EnviroShape=0;
SharingNess=[0,0.1,10^-3];
ChangeTime=50;
baseFit=1;
baseDeath=0.35;
FitDeltaLin=0.0;
FitDeltaDeath=baseDeath/5;
FitDeltaCombo=0.3;
fileName='GrandSet4.mat';

InitialRatios=[0.54,0.06]


Jumps=[1,-1,2^Q,-2^Q];

numRepeats=100;
results=zeros(20,2);
timeTime=zeros(20,1);
Timecourse=zeros(20,1);

AverageTimecourse=zeros(1,20);
AverageWinTimecourse=zeros(1,20);
AverageLoseTimecourse=zeros(1,20);

AverageFinalDense=zeros(1,20);
AverageWinFinalDense=zeros(1,20);
AverageLoseFinalDense=zeros(1,20);

percentages=zeros(3,20);

comments={};
%load(fileName)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q=1;
baseFit=1;
baseDeath=0.35;
FitDeltaLin=0.0;
FitDeltaCombo=0.3;

SimCondition=1;
comments{1}= 'baseline';

for(iii=(1:numRepeats))
    tic()
 [Trajectory] = PlasmidSpreadFunction2d_BinaryTree(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= size(Trajectory,2);

results(iii,:)=Trajectory(:,end)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0),2);
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0),1);

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SimCondition=15;
comments{15}= 'Diagonal Jump';
Jumps=[1,-1,2^Q,-2^Q,2^Q+1,-2^Q+1,2^Q-1,-2^Q-1];

for(iii=(1:numRepeats))
    tic()
[Trajectory] = PlasmidSpreadFunction2d_BinaryTree(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= size(Trajectory,2);

results(iii,:)=Trajectory(:,end)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0),2);
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0),1);

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);

Jumps=[1,-1,2^Q,-2^Q];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SimCondition=2;

q=3;
baseFit=1;
FitDeltaLin=0.1;
FitDeltaCombo=0.0;
mutationRate=10^-4;
Q=7;

for(iii=(1:numRepeats))
    tic()
[Trajectory] = PlasmidSpreadFunction2d_BinaryTree(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= size(Trajectory,2);

results(iii,:)=Trajectory(:,end)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0),1);
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0),1);

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);

comments{2}= '3 Plasmid; linear sum';

%%%%%%%%%%%%%%%%%%%%%%5
SimCondition=3;

q=3;
baseFit=1;
FitDeltaLin=0.0;
FitDeltaCombo=0.3;
mutationRate=10^-4;
Q=7;

for(iii=(1:numRepeats))
    tic()
[Trajectory] = PlasmidSpreadFunction2d_BinaryTree(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= size(Trajectory,2);

results(iii,:)=Trajectory(:,end)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0),1);
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0),1);

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);

comments{3}= '3 Plasmid; No correlation';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SimCondition=4;

q=3;
baseFit=1;
FitDeltaLin=0.05;
FitDeltaCombo=0.15;
mutationRate=10^-4;
Q=7;

for(iii=(1:numRepeats))
    tic()
[Trajectory] = PlasmidSpreadFunction2d_BinaryTree(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= size(Trajectory,2);

results(iii,:)=Trajectory(:,end)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0),1);
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0),1);

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);

comments{4}= '3 Plasmid; Partial Correlation';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

comments{5}= '1 Plasmid; Block Environment';
SimCondition=5;

q=1;
baseFit=1;
FitDeltaLin=0.0;
FitDeltaCombo=0.3;
EnviroShape=1; %%Two block environment

for(iii=(1:numRepeats))
    tic()
[Trajectory] = PlasmidSpreadFunction2d_BinaryTree(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= size(Trajectory,2);

results(iii,:)=Trajectory(:,end)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));
percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% 

comments{6}= '1 Plasmid; Checkerboard Environment';
SimCondition=6;

q=1;
baseFit=1;
FitDeltaLin=0.0;
FitDeltaCombo=0.3;
EnviroShape=2; %%Checkboard

for(iii=(1:numRepeats))
    tic()
[Trajectory] = PlasmidSpreadFunction2d_BinaryTree(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= size(Trajectory,2);

results(iii,:)=Trajectory(:,end)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0),2);
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0),1);

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

comments{7}= '1 Plasmid; Random Environment';
SimCondition=7;

q=1;
baseFit=1;
FitDeltaLin=0.0;
FitDeltaCombo=0.3;
EnviroShape=3; %%random environment

for(iii=(1:numRepeats))
    tic()
[Trajectory] = PlasmidSpreadFunction2d_BinaryTree(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= size(Trajectory,2);

results(iii,:)=Trajectory(:,end)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0),2);
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0),1);

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

comments{8}= '1 Plasmid; Fast turnover (T=5)';
SimCondition=8;

q=1;
baseFit=1;
FitDeltaLin=0.0;
FitDeltaCombo=0.3;
EnviroShape=0; %%No Environment
ChangeTime=5;

for(iii=(1:numRepeats))
    tic()
[Trajectory] = PlasmidSpreadFunction2d_BinaryTree(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= size(Trajectory,2);

results(iii,:)=Trajectory(:,end)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0),2);
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0),1);

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
comments{9}= '1 Plasmid; Slow turnover (T=250)';
SimCondition=9;

q=1;
baseFit=1;
FitDeltaLin=0.0;
FitDeltaCombo=0.3;
EnviroShape=0; %%No Environment
ChangeTime=250;

for(iii=(1:numRepeats))
    tic()
[Trajectory] = PlasmidSpreadFunction2d_BinaryTree(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= size(Trajectory,2);

results(iii,:)=Trajectory(:,end)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0),2);
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0),1);

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

comments{10}= '1 Plasmid; Weak Restriction (10^-2)';
SimCondition=10;

q=1;
baseFit=1;
FitDeltaLin=0.0;
FitDeltaCombo=0.3;
EnviroShape=0; %%No Environment
ChangeTime=50;
SharingNess=[0,0.1,10^-2];
for(iii=(1:numRepeats))
    tic()
[Trajectory] = PlasmidSpreadFunction2d_BinaryTree(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= size(Trajectory,2);

results(iii,:)=Trajectory(:,end)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0),2);
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0),1);

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

comments{11}= '1 Plasmid; Harsh Restriction (10^-5)';
SimCondition=11;

q=1;
baseFit=1;
FitDeltaLin=0.0;
FitDeltaCombo=0.3;
EnviroShape=0; %%No Environment
ChangeTime=50;
SharingNess=[0,0.1,10^-5];

for(iii=(1:numRepeats))
    tic()
[Trajectory] = PlasmidSpreadFunction2d_BinaryTree(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= size(Trajectory,2);

results(iii,:)=Trajectory(:,end)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0),2);
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0),1);

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


comments{14}= 'Control: c(S)=c(G)';
SimCondition=14;

q=1;
baseFit=1;
FitDeltaLin=0.0;
FitDeltaCombo=0.3;
EnviroShape=0; %%No Environment
ChangeTime=50;
SharingNess=[0,0.1,0.1];

for(iii=(1:numRepeats))
    tic()
[Trajectory] = PlasmidSpreadFunction2d_BinaryTree(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= size(Trajectory,2);

results(iii,:)=Trajectory(:,end)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0),2);
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0),1);

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

comments{12}= '1 Plasmid; Higher deathrate, lower density (base death=0.75)';
SimCondition=12;

q=1;
baseFit=1;
FitDeltaLin=0.0;
FitDeltaCombo=0.3;
EnviroShape=0; %%No Environment
ChangeTime=50;
baseDeath=0.60;
FitDeltaDeath=baseDeath/5;

for(iii=(1:numRepeats))
    tic()
[Trajectory] = PlasmidSpreadFunction2d_BinaryTree(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= size(Trajectory,2);

results(iii,:)=Trajectory(:,end)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0),2);
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0),1);

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

comments{20}= '1 Plasmid; Lower deathrate, higher density (base death=0.05) ACCELERATED';
SimCondition=20;

q=1;
baseFit=1;
FitDeltaLin=0.0;
FitDeltaCombo=0.3;
EnviroShape=0; %%No Environment
ChangeTime=50;
baseDeath=0.05;
FitDeltaDeath=baseDeath/5;
SharingNess=[0,0.1,10^-3];

for(iii=(1:numRepeats))
    tic()
[Trajectory] = PlasmidSpreadFunction2d_BinaryTreePheonix(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= Trajectory(3);

results(iii,:)=Trajectory(1:2)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0,2));
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0,1));

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);
save(fileName); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


comments{16}= 'Environment; Gradient';
SimCondition=16;

SharingNess=[0,0.1,10^-3];
q=1;
k=4;
baseFit=1;
FitDeltaLin=0.0;
FitDeltaCombo=0.3;
EnviroShape=2; %%Gradient
ChangeTime=50;
baseDeath=0.35;
FitDeltaDeath=baseDeath/5;

for(iii=(1:numRepeats))
    tic()
[Trajectory] = PlasmidSpreadFunction2d_BinaryTreeGradient(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= size(Trajectory,2);

results(iii,:)=Trajectory(:,end)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0,2));
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0,1));

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);
save(fileName); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q=1;
baseFit=1;
baseDeath=0.35;
FitDeltaLin=0.0;
FitDeltaCombo=0.3;
SharingNess=[0,0.1,10^-3];

SimCondition=17;
comments{17}= 'High Mutation (x30)';

for(iii=(1:numRepeats))
    tic()
[Trajectory] = PlasmidSpreadFunction2d_BinaryTree(Q,k,q,mutationRate*30,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= size(Trajectory,2);

results(iii,:)=Trajectory(:,end)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0,2));
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0,1));

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);
save(fileName); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q=1;
baseFit=1;
baseDeath=0.35;
FitDeltaLin=0.0;
FitDeltaCombo=0.3;

SimCondition=18;
comments{18}= 'Low Mutation Rate';
SharingNess=[0,0.1,10^-3];

for(iii=(1:numRepeats))
    tic()
[Trajectory] = PlasmidSpreadFunction2d_BinaryTree(Q,k,q,mutationRate/30,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= size(Trajectory,2);

results(iii,:)=Trajectory(:,end)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0,2));
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0,1));

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);
save(fileName); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q=3;
k=25;
baseFit=1;
baseDeath=0.35;
FitDeltaLin=0.05;
FitDeltaCombo=0.15;

SimCondition=19;
comments{19}= 'q=3, k=25';
SharingNess=[0,0.1,10^-3];

for(iii=(1:numRepeats))
    tic()
[Trajectory] = PlasmidSpreadFunction2d_BinaryTreeNoRecord(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= Trajectory(3);

results(iii,:)=Trajectory(1:2)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0,2));
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0,1));

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save(fileName);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
q=1;
k=4;
baseFit=1;
baseDeath=0.35;
FitDeltaLin=0.05;
FitDeltaCombo=0.15;

SimCondition=21;
comments{21}= 'N=128';
SharingNess=[0,0.1,10^-3];
Q=7; %N=128;

for(iii=(1:numRepeats))
    tic()
[Trajectory] = PlasmidSpreadFunction2d_BinaryTreeNoRecord(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= Trajectory(3);

results(iii,:)=Trajectory(1:2)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0,2));
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0,1));

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save(fileName);

q=1;
k=4;
baseFit=1;
baseDeath=0.35;
FitDeltaLin=0.05;
FitDeltaCombo=0.15;

SimCondition=22;
comments{22}= 'Plasmid PULL simulation';
SharingNess=[0,0.1,10^-3];
Q=6; %N=64;

for(iii=(1:numRepeats))
    tic()
[Trajectory] = PlasmidSpreadFunction2d_BinaryTree_Pull(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= size(Trajectory,2);
results(iii,:)=Trajectory(:,end)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0,2));
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0,1));

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save(fileName);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q=1;
baseFit=1;
baseDeath=0.35;
FitDeltaLin=0.0;
FitDeltaCombo=0.3;

SimCondition=23;
comments{23}= 'GroupSelection; regularTiming';

for(iii=(1:numRepeats))
    tic()
 [Trajectory] = PlasmidSpreadFunction2d_BinaryTree_GroupSelection(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= size(Trajectory,2);

results(iii,:)=Trajectory(:,end)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0),2);
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0),1);

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);

q=1;
baseFit=1;
baseDeath=0.35;
FitDeltaLin=0.0;
FitDeltaCombo=0.3;
numRepeats=400;

SimCondition=24;
comments{24}= 'GroupSelection; MORE REPEATS';

for(iii=(1:numRepeats))
    tic()
 [Trajectory] = PlasmidSpreadFunction2d_BinaryTree_GroupSelection(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= size(Trajectory,2);

results(iii,:)=Trajectory(:,end)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0),1);
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0),1);

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q=1;
baseFit=1;
baseDeath=0.35;
FitDeltaLin=0.0;
FitDeltaCombo=0.3;
numRepeats=100;

SimCondition=25;
comments{25}= 'GroupSelection; Traulsen implementation';

for(iii=(1:numRepeats))
    tic()
 [Trajectory] = PlasmidSpreadFunction2d_GroupSelection_2nd_implementation(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= size(Trajectory,2);

results(iii,:)=Trajectory(:,end)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0),1);
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0),1);

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


q=1;
baseFit=1;
baseDeath=0.35;
FitDeltaLin=0.0;
FitDeltaCombo=0.3;
numRepeats=100;

SimCondition=25;
comments{25}= 'GroupSelection; Traulsen implementation, no Time Limit, 10% initial';

InitialRatios=[0.54,0.06]


for(iii=(1:numRepeats))
    tic()
 [Trajectory] = PlasmidSpreadFunction2d_GroupSelection_Traulsen_noRecord(Q,k,q,mutationRate,Jumps,EnviroShape,SharingNess,baseFit,baseDeath,FitDeltaLin,FitDeltaDeath,FitDeltaCombo,1000*SimCondition+iii,ChangeTime,InitialRatios);

timeTime(iii)=toc();
Timecourse(iii)= Trajectory(3);

results(iii,:)=Trajectory(1:2)';
save(fileName);

end

AverageTimecourse(SimCondition)= mean(Timecourse);
AverageWinTimecourse(SimCondition)= mean(Timecourse(results(:,1)==0));
AverageLoseTimecourse(SimCondition)= mean(Timecourse(results(:,2)==0));

AverageFinalDense(SimCondition)= mean(results(:,2)+results(:,1));
AverageWinFinalDense(SimCondition)= mean(results(results(:,1)==0),1);
AverageLoseFinalDense(SimCondition)= mean(results(results(:,2)==0),1);

percentages(1,SimCondition)= mean(results(:,1)==0);
percentages(2,SimCondition)= mean(results(:,2)==0);
percentages(3,SimCondition)= 1-percentages(1,SimCondition)-percentages(2,SimCondition);
