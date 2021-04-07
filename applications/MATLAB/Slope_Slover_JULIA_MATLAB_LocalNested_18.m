function [ ] = Slope_Slover_JULIA_MATLAB_LocalNested(ProcID)
% ProcID=1
% NumberOfGaussPoints=[7 16 28 37];



%%%%%%%%%SELECT CORRECT SYSTEM
str_folder="/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/dev/MultilevelEstimators/applications/SPDE/data/";
%str_folder="/home/philippeb/.julia/packages/MultilevelEstimators/l8j9n/applications/SPDE/data/";
%str_folder="/home/philippe/.julia/dev/Applications/applications/SPDE/data/";
%%%%%%%%%%%%%



str_interm=str_folder+"Interm28/Slope/";


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
%Gausspoints{1}=[6,8,11,13,3,1,4,7,2,9,10,12,5];%L0_M0_S0

Gausspoints{1}=[7,16,10,15,2,11,3,1,9,6,13,4,12,8,14,5];%L1_M1_S1
%Gausspoints{3}=[6,8,3,11,13,1,4,7,9,10,2,12,5];%L0_M1_S0 %done

Gausspoints{2}=[13,17,19,10,5,16,6,1,9,12,18,7,8,14,15,11,2,3,4];%L2_M2_S2
Gausspoints{3}=[7,16,15,10,2,11,3,1,9,6,13,4,8,12,14,5];%L1_M2_S1 %done

Gausspoints{4}=[4,20,8,22,13,14,19,9,15,1,12,3,21,16,11,17,10,18,2,26,28,7,25,6,27,5,23,24];%L3_M3_S3 %done
Gausspoints{5}=[13,17,2,19,10,5,16,3,6,1,9,12,18,7,8,14,4,15,11];%L2_M3_S2 %done


Gausspoints{6}=[4,21,24,27,13,14,30,10,26,5,20,15,6,1,9,12,3,25,7,31,8,28,16,11,29,22,23,2,33,17,36,32,18,37,34,19,35];%L4_M4_S4 %done
Gausspoints{7}=[4,26,28,20,7,8,22,13,19,14,25,9,15,1,12,6,3,27,16,21,11,17,10,5,18,23,24,2];%L3_M4_S3 %done

Gausspoints{8}=[9,7,18,2,25,12,1,10,3,37,34,11,29,28,6,46,43,38,35,30,8,17,24,39,44,33,50,52,36,42,45,49,55,51,59,60,61,4,16,13,27,22,20,19,31,5,21,14,23,15,40,32,26,47,41,48,53,57,54,58,56]; %L5_M5_S5 %done
Gausspoints{9}=[4,21,24,27,13,33,14,17,26,30,36,32,10,5,20,15,18,6,1,9,3,25,12,7,8,37,34,28,31,19,16,35,11,29,22,23,2];%L4_M5_S4 %done


Gausspoints{10}=[8,6,19,1,17,30,20,3,14,35,11,31,26,4,22,15,41,44,37,46,33,25,2,21,5,52,50,47,42,36,32,16,7,18,34,29,40,48,59,51,57,54,43,49,45,55,64,60,53,62,56,58,61,69,63,66,70,68,71,72,73,12,9,28,23,39,10,13,24,27,38,65,67]; %L6_M6_S6 %done
Gausspoints{11}=[9,7,18,4,16,27,13,2,12,25,10,22,20,3,19,11,31,34,29,37,28,21,5,14,6,43,40,38,35,30,23,15,8,17,24,26,32,39,52,44,50,47,33,41,36,48,46,53,42,55,49,45,51,57,1,54,58,56,59,60,61]; %L5_M6_S5 %done


Gausspoints{12}=[29,49,61,52,82,85,67,64,20,88,103,73,17,106,70,1,97,109,102,78,4,14,76,112,72,114,7,66,40,108,90,54,3,6,9,15,18,21,30,80,50,86,10,104,39,62,8,110,68,111,75,13,5,74,38,113,2,69,105,16,41,77,107,87,19,63,71,89,81,51,65,53,28,55,58,26,79,31,23,43,91,37,34,115,42,118,94,46,48,121,36,96,124,126,100,11,120,84,60,12,24,27,56,116,99,122,123,92,32,125,44,98,45,93,117,101,119,33,95,47,22,35,83,25,57,59]; %L7_M7_S7
Gausspoints{13}=[8,6,1,19,17,12,9,30,23,28,20,14,31,35,39,11,3,26,4,15,22,37,44,41,10,25,33,2,46,21,13,5,52,50,47,36,32,24,7,16,18,27,42,34,59,29,48,40,38,54,57,51,49,43,45,55,53,64,62,60,58,56,61,67,65,69,63,66,70,71,68,72,73]; %L6_M7_S6 5done

%Gausspoints{14}=[36,56,59,68,89,110,50,71,113,74,27,140,98,62,53,143,24,152,38,92,146,77,134,155,101,2,21,164,158,43,65,97,149,167,139,18,67,170,8,163,79,103,151,137,95,161,41,173,55,15,11,175,169,73,157,145,115,61,4,10,13,16,19,22,25,28,37,57,87,111,141,153,69,165,171,12,14,51,172,40,160,94,136,147,99,75,159,9,174,63,17,135,166,148,93,64,39,162,168,20,3,100,154,138,76,150,96,42,156,23,142,52,66,102,144,26,78,112,70,54,114,88,72,58,60,35,44,33,86,47,30,104,128,107,80,116,131,122,127,5,119,85,83,121,133,109,49,125,91,1,7,31,34,124,45,105,129,117,82,81,118,6,123,126,130,120,84,106,132,108,29,46,90,32,48]; %L8_M8_S8
%Gausspoints{15}=[29,55,58,49,26,79,61,52,82,31,23,85,67,43,64,88,20,103,37,73,91,34,97,106,70,1,17,115,109,78,46,42,94,118,102,14,48,121,4,114,36,72,96,100,76,112,40,124,66,11,7,126,120,54,108,90,84,60,3,6,9,12,15,18,21,24,30,56,27,80,86,104,50,116,122,8,10,62,123,75,111,39,99,92,68,32,110,5,125,44,13,98,117,93,74,45,38,113,119,16,2,69,105,101,33,95,77,41,107,19,87,63,47,71,89,22,35,81,51,65,83,25,53,57,59,28];%L7_M8_S7

%elementPlane={'plane3','plane6','plane10','plane15','plane21','plane28','plane36t','plane45','plane55'};
elementPlane={'plane6','plane10','plane15','plane21','plane28','plane36t','plane45'};
keySet={'M0_S0','M1_S1','M1_S0','M2_S2','M2_S1','M3_S3','M3_S2','M4_S4','M4_S3','M5_S5','M5_S4','M6_S6','M6_S5'};
OrderDic=[1,2,1,3,2,4,3,5,4,6,5,7,6];
val=[1,2,3,4,5,6,7];
%%

% NumberOfGaussPoints=[6 7 13 16 19 28 37 61 73]; %Paper
NumberOfGaussPoints=[16 19 28 37 61 73 126]; %Paper
%7
%7 16
%16 19
%19 28
%28 37
%37 61
%61 73
%61 73
%61 73

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

[elempos,col]=find(Elements==CornerNode(1));
elem=max(elempos);

NodesElem=Elements(elem,5:end);

% posNodes=Nodes(NodesElem,2:3);
% 
% xcoord=sum(posNodes(:,1))/3;
% ycoord=sum(posNodes(:,2))/3;


posNodesTriangleLocal=[[0 1 0]',[0 0 1]'];
xcoord=sum(posNodesTriangleLocal(:,1))/3;
ycoord=sum(posNodesTriangleLocal(:,2))/3;

 Solcenter=-sh_t6(xcoord,ycoord)'*selectdof(DOF,NodesElem+0.02)*U(:,end); 



% res1=-selectdof(DOF,CornerNode(1)+0.02)*[zeros(size(DOF)) U];
%res2=-selectdof(DOF,9+0.02)*[zeros(size(DOF)) U];
%res=[res1 res2];
% res=res1(end);

dlmwrite(fileRES,Solcenter,'delimiter',' ','precision',15);


end
