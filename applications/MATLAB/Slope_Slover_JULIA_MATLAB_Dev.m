function [ ] = Slope_Slover_JULIA_MATLAB_Dev(ProcID)
% ProcID=1
% NumberOfGaussPoints=[7 16 28 37]; 
NumberOfGaussPoints=[7 13  19 25 28 33 37 61 73]; %Paper


%%%%%%%%%SELECT CORRECT SYSTEM
str_folder="/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/packages/MultilevelEstimators/KoRDo/applications/SPDE/data/";
%str_folder="/home/philippeb/.julia/packages/MultilevelEstimators/l8j9n/applications/SPDE/data/";
%str_folder="/home/philippe/.julia/dev/Applications/applications/SPDE/data/";
%%%%%%%%%%%%%



str_interm=str_folder+"Interm/Slope/";


fileC=strcat(str_interm,strcat("Cohesion_",num2str(ProcID)));
fileRES=strcat(str_interm,strcat("Res_",num2str(ProcID)));
fileLevel=strcat(str_interm,strcat("Index_",num2str(ProcID)));
filedoRestriction=strcat(str_interm,strcat("doRestriction_",num2str(ProcID)));
dores=dlmread(filedoRestriction);



level=dlmread(fileLevel);
level_num=level;
HigherOrderHandle=strcat(str_interm,strcat("HighOrder_",num2str(ProcID)));
HigherOrder=dlmread(HigherOrderHandle);

if(length(level)==1)
    
    str=strcat('_',num2str(level(1)));
    level=str;
    
    
    if(HigherOrder==0) %Pure h refinement
        Order=1;
        
        str_elements=str_folder+"Mesh/Slope/h_refinement/";
        
        
        files_elements=strcat(str_elements,"Elements_L",level,".txt");
        Elements=dlmread(files_elements);
        
        nElem=size(Elements,1);
    elseif(HigherOrder==1) %Pure p refinement
        Order=level_num+1;
        str_elements=str_folder+"Mesh/Slope/p_refinement/";
        
        
        files_elements=strcat(str_elements,"Elements_L",level,".txt");
        Elements=dlmread(files_elements);
        
        nElem=size(Elements,1);
        EltOpts.GP_k_dp=1;
        %nelem TODO
    end
else
    
    if(HigherOrder==3) %MultiIndex hp refinement
        Order=level(2)+1;
        str_elements=str_folder+"Mesh/hp_refinement/Slope/";
        str='_';
        for len=1:length(level)
            str=strcat(str,num2str(level(len)));
            
        end
        
        files_elements=strcat(str_elements,"Elements_L",str,".txt");
        Elements=dlmread(files_elements);
        
        nElem=size(Elements,1);
        EltOpts.GP_k_dp=1;
        level=str;
        
    end
end


fileExt=strcat(str_interm,"Extension","_",num2str(ProcID));
fx=fopen(fileExt);
Extension=fscanf(fx, '%s');
fclose(fx);
Extension=string(Extension);

Gausspoints={};

%% For specific runs
% Gausspoints{1}=[5,3,1,4,7,2,6];
% Gausspoints{2}=[7,2,3,1,6,4,5,16,10,15,11,9,13,12,8,14];
% Gausspoints{3}=[5,3,4,1,7,2,6];
% Gausspoints{4}=[4,20,22,13,14,19,15,1,12,3,21,16,11,17,18,2,26,28,7,8,25,9,6,27,10,5,23,24];
% Gausspoints{5}=[7,16,15,10,2,11,3,1,9,6,13,4,8,12,14,5];
% Gausspoints{6}=[4,21,24,27,13,14,30,10,26,5,20,15,6,1,9,12,3,25,7,31,8,28,16,11,29,22,23,2,33,17,36,32,18,37,34,19,35];
% Gausspoints{7}=[4,26,28,20,7,8,22,13,19,14,25,9,15,1,12,6,3,27,16,21,11,17,10,5,18,23,24,2];
% elementPlane={'plane3','plane10','plane15','plane21'};
% keySet={'M0_S0','M1_S1','M1_S0','M2_S2','M2_S1','M3_S3','M3_S2'};
% OrderDic=[1,2,1,3,2,4,3];
% val=[1,2,3,4];
Gausspoints{1}=[5,3,1,4,7,2,6];

Gausspoints{2}=[6,8,11,1,7,10,12,13,3,4,2,9,5];
Gausspoints{3}=[5,3,4,1,7,6,2];

Gausspoints{4}=[13,17,19,10,16,1,9,12,18,8,14,15,11,2,5,3,6,7,4];
Gausspoints{5}=[6,8,11,3,13,1,9,7,4,12,10,2,5];

Gausspoints{6}=[2,3,13,5,14,4,11,18,16,9,6,15,12,19,20,17,22,21,24,7,8,10,1,23,25];%L3_M3_S3 %done
Gausspoints{7}=[13,17,19,2,10,16,5,3,1,9,12,7,18,8,6,4,14,15,11];%L2_M3_S2 %done


Gausspoints{8}=[4,28,7,20,8,22,13,19,9,15,1,12,6,3,27,21,16,11,17,10,18,5,23,24,2,26,14,25];%L4_M4_S4 %done
Gausspoints{9}=[2,7,8,3,5,13,10,4,18,14,16,11,9,1,6,12,15,20,22,17,21,23,24,19,25];%L3_M4_S3 %done

Gausspoints{10}=[14,28,31,22,11,16,1,25,4,27,8,21,33,3,6,12,15,29,23,5,24,2,20,10,26,30,32,13,19,9,17,7,18]; %L5_M5_S5 %done
Gausspoints{11}=[4,26,28,20,7,13,8,22,14,19,1,12,25,9,15,6,3,27,21,16,17,10,11,5,18,23,24,2];%L4_M5_S4 %done


Gausspoints{12}=[4,21,24,27,13,33,14,30,36,10,26,5,32,20,15,6,9,12,3,25,37,7,31,8,34,28,16,35,11,29,22,23,2,17,1,18,19]; %L6_M6_S6 %done
Gausspoints{13}=[14,28,31,22,11,16,1,25,19,8,27,4,21,33,3,6,9,12,15,29,17,5,23,7,18,24,2,20,10,26,30,32,13]; %L5_M6_S5 %done


Gausspoints{14}=[9,7,18,2,25,12,1,10,3,37,34,11,29,28,6,46,43,38,35,30,8,17,24,39,44,33,50,52,36,42,45,49,55,51,59,60,61,4,16,13,27,22,20,19,31,5,21,14,23,15,40,32,26,47,41,48,53,57,54,58,56]; %L7_M7_S7
Gausspoints{15}=[4,21,24,27,13,33,14,17,26,30,36,32,10,5,20,15,18,6,1,9,3,25,12,7,8,37,34,28,31,19,16,35,11,29,22,23,2]; %L6_M7_S6 5done

Gausspoints{16}=[8,6,19,1,17,30,20,3,14,35,11,31,26,4,22,15,41,44,37,46,33,25,2,21,5,52,50,47,42,36,32,16,7,18,34,29,40,48,59,51,57,54,43,49,45,55,64,60,53,62,56,58,61,69,63,66,70,68,71,72,73,12,9,28,23,39,10,13,24,27,38,65,67]; %L8_M8_S8
Gausspoints{17}=[9,7,18,4,16,27,13,2,12,25,10,22,20,3,19,11,31,34,29,37,28,21,5,14,6,43,40,38,35,30,23,15,8,17,24,26,32,39,52,44,50,47,33,41,36,48,46,53,42,55,49,45,51,57,1,54,58,56,59,60,61];%L7_M8_S7

elementPlane={'plane3','plane6','plane10','plane15','plane21','plane28','plane36t','plane45','plane55'};
keySet={'M0_S0','M1_S1','M1_S0','M2_S2','M2_S1','M3_S3','M3_S2','M4_S4','M4_S3','M5_S5','M5_S4','M6_S6','M6_S5','M7_S7','M7_S6','M8_S8','M8_S7'};
OrderDic=[1,2,1,3,2,4,3,5,4,6,5,7,6,8,7,9,8];
val=[1,2,3,4,5,6,7,8,9];
%%




Dict1 = containers.Map(keySet,OrderDic);
DicPlane=containers.Map(val,elementPlane);

DicGaussPoints=containers.Map(keySet,Gausspoints);





fileC=strcat(str_interm,"Cohesion",level,Extension,"_",num2str(ProcID));
file_nodes=strcat(str_elements,"Nodes_L",level,".txt");



c=dlmread(fileC);

if(HigherOrder==0)
    Nodes=dlmread(file_nodes);
    file_nodes_centers=strcat(str_elements,"ElementsCenter_L",level,".txt");
    nodes_centers=dlmread(file_nodes_centers);
  %  c=c(1:nElem);
    c=griddata(Nodes(:,2),Nodes(:,3),c,nodes_centers(:,1),nodes_centers(:,2),'cubic');  % Cohesion in N/m^2
% F=scatteredInterpolant(Nodes(:,2),Nodes(:,3),c);
% c=F(nodes_centers(:,1),nodes_centers(:,2));

else
    
    %c=c(1:NumberOfGaussPoints(Order)*nElem);
    %  c=reshape(c,nElem,NumberOfGaussPoints(Order));
    c=transpose(reshape(c,NumberOfGaussPoints(Dict1(char(Extension))),nElem));
end




file_nodes=strcat(str_elements,"Nodes_L",level,".txt");
files_elements=strcat(str_elements,"Elements_L",level,".txt");
Elements=dlmread(files_elements);
Nodes=dlmread(file_nodes);



%%%%%%%%%%

%% Parameters
E=30e6;         % Young's modulus in N/m^2
nu=0.25;        % Poisson ratio
rho=1330;       % Density [kg?m^3]
phi=20;         % Friction angle in degrees

% Match Mohr-Coulomb at shear meridian
alpha_dp=sind(phi)/3;
k_dp=c.*cosd(phi);

%  Extension="M33";
EltOpts.nXi=NumberOfGaussPoints(Dict1(char(Extension)));
elementPlane=DicPlane(Dict1(char(Extension)));

EltOpts.Sort_nXi=DicGaussPoints(char(Extension));

% switch Order
%     case 1
%         elementPlane='plane3';
%         EltOpts.nXi=NumberOfGaussPoints(1);
%         EltOpts.Sort_nXi=[5,1,3,4,2,7,6];
%     case 2
%         elementPlane='plane10';
%         EltOpts.nXi=NumberOfGaussPoints(2);
%          EltOpts.Sort_nXi=;
% 
%     case 3
%         elementPlane='plane21';
%         EltOpts.nXi=NumberOfGaussPoints(3);
%         EltOpts.Sort_nXi=;
% 
%     case 4
%         elementPlane='plane36t';
%         EltOpts.nXi=NumberOfGaussPoints(4);
%         EltOpts.Sort_nXi=;
%  
% end



StepSize=0.1/10;


%% Create model
EltOpts.problem='2dstrain';
EltOpts.nl='elastoplastic';
Types={1 elementPlane EltOpts};
Materials = cell(nElem,3);



nElem=size(Elements,1);



 for k =1:nElem, Materials(k,:)={k 'elastoplastic' {'isotropic' [E nu rho] 'dp' [alpha_dp k_dp(k,:)] 'multilinear' [1 k_dp(k,:)+E/300]}};end

Sections=[1];
%




LeftNodes=selectnode(Nodes,-1e-6,-inf,-inf,1e-6,inf,inf);
RightNodes=selectnode(Nodes,max(Nodes(:,2))-1e-6,-inf,-inf,max(Nodes(:,2))+1e-6,inf,inf);
BottomNodes=selectnode(Nodes,-inf,-1e-6,-inf,inf,1e-6,inf);



% BottomNodes2=model.getNodes('lines',1);
% TopNodes=model.getNodes('lines',3);

% Create list of DOFs
DOF=getdof(Elements,Types);
seldof=[LeftNodes(:,1)+0.01;RightNodes(:,1)+0.01;BottomNodes(:,1)+0.01;BottomNodes(:,1)+0.02];
InitalDof=DOF;
DOF=removedof(DOF,seldof);


% Load vector
DLoads=accel([0 9.81 0],Elements,Types,Sections,{1 E nu rho});
P=elemloads(DLoads,Nodes,Elements,Types,DOF);

%% Compute equilibrium path



CornerNode=selectnode(Nodes,15,14,0); % QOI



Options=struct;

Options.verbose=false;
Div=1;
% while (Div==1)
    Pload=(0:StepSize:1); %load
     [U,~,~,~,Div]=solver_nr(Nodes,Elements,Types,Sections,Materials,DOF,P,Pload,Options,[],[],[],[],[]);
%     
%     res=-selectdof(DOF,8+0.02)*[zeros(size(DOF)) U];
%     res=res(end);
%     count=1;
%     if(Div==1||res>1||res<0)
%         StepSize=StepSize/2;
%         Div=1;
%         count=count+1;
%         if(count>10)
%             
%             Div=0;
%             U=zeros(length(Nodes),1);
%         end
%     end
%     
%     
% end
res1=-selectdof(DOF,CornerNode(1)+0.02)*[zeros(size(DOF)) U];
%res2=-selectdof(DOF,9+0.02)*[zeros(size(DOF)) U];
%res=[res1 res2];
res=res1(end);

dlmwrite(fileRES,res,'delimiter',' ','precision',15);


end

