classdef VREMF < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained>
% Coevolutionary constrained multi-objective optimization framework
% type --- 1 --- Type of operator (1. GA 2. DE)

%------------------------------- Reference --------------------------------
% Y. Tian, T. Zhang, J. Xiao, X. Zhang, and Y. Jin, A coevolutionary
% framework for constrained multi-objective optimization problems, IEEE
% Transactions on Evolutionary Computation, 2021, 25(1): 102-116.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population1   = Problem.Initialization();
            Population2   = Problem.Initialization();
            Fitness1      = CalFitness(Population1.objs,Population1.cons);
            Fitness2      = CalFitness(Population2.objs);
            threshold1    = 1e-6;                           
            Gen           = 2;                              
            Objvalues(1)  = sum(sum(Population2.objs,1));   
            Objvalues(2)  = 1000;
            result        = 0;                              
            Con           = [];
            Dec           = [];
            T             = 100;
            %% Optimization
            while Algorithm.NotTerminated(Population1)
               if  result == 1                  
                   if  rf < 1
                       [OS, UOS] = CACK(Problem, Con, Dec);
                       B = [min(Population2.decs);max(Population2.decs);min(Population1.decs);max(Population1.decs)];  
                       MatingPool1 = TournamentSelection(2,2*Problem.N,Fitness1);
                       MatingPool2 = TournamentSelection(2,2*Problem.N,Fitness2);
                       Offspring1  = OperatorDE(Problem,Population1,Population1(MatingPool1(1:end/2)),Population1(MatingPool1(end/2+1:end)));
                       Offspring2  = OperatorDE2(Problem,Population2,Population2(MatingPool2(1:end/2)),Population2(MatingPool2(end/2+1:end)),B,OS,UOS);
                       [Population1,Fitness1] = EnvironmentalSelection([Population1,Offspring1,Offspring2],Problem.N,true);
                       [Population2,Fitness2] = EnvironmentalSelection([Population2,Offspring1,Offspring2],Problem.N,true);
                   else
                       MatingPool1 = TournamentSelection(2,2*Problem.N,Fitness1);
                       MatingPool2 = TournamentSelection(2,2*Problem.N,Fitness2);
                       Offspring1  = OperatorDE(Problem,Population1,Population1(MatingPool1(1:end/2)),Population1(MatingPool1(end/2+1:end)));
                       Offspring2  = OperatorDE1(Problem,Population2,Population2(MatingPool2(1:end/2)),Population2(MatingPool2(end/2+1:end)),Population1.decs);
                       [Population1,Fitness1] = EnvironmentalSelection([Population1,Offspring1,Offspring2],Problem.N,true);
                       [Population2,Fitness2] = EnvironmentalSelection([Population2,Offspring1,Offspring2],Problem.N,false);
                   end
               else                        
                   MatingPool1 = TournamentSelection(2,2*Problem.N,Fitness1);
                   MatingPool2 = TournamentSelection(2,2*Problem.N,Fitness2);
                   Offspring1  = OperatorDE(Problem,Population1,Population1(MatingPool1(1:end/2)),Population1(MatingPool1(end/2+1:end)));
                   Offspring2  = OperatorDE1(Problem,Population2,Population2(MatingPool2(1:end/2)),Population2(MatingPool2(end/2+1:end)),Population1.decs);
                   [Population1,Fitness1] = EnvironmentalSelection([Population1,Offspring1,Offspring2],Problem.N,true);
                   [Population2,Fitness2] = EnvironmentalSelection([Population2,Offspring1,Offspring2],Problem.N,false);
               end
               if result == 0
                   result          = is_update(Objvalues,Gen,Population2,Problem.N,threshold1,Problem.M);
               end
               Gen               = Gen+1;
               Objvalues(Gen)    = sum(sum(abs(Population2.objs),1));                                      
               CV        = sum(max(0,Population1.cons),2);
               rf        = sum(CV <= 1e-6) / length(Population1);                                         
               [m, ~]    = size(Con);
               if m > T
                   Con(m-T,:)   = [];
                   Dec(m-T,:)   = [];
                   Con          = [Con;mean(Population2.cons)];
                   Dec          = [Dec;mean(Population2.decs)];
               else
                   Con          = [Con;mean(Population2.cons)];
                   Dec          = [Dec;mean(Population2.decs)];
               end
            end
        end
    end
end