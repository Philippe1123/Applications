function [] = Populate_3d_MATLAB( maxlevel )
%POPULATE_3D_MATLAB Summary of this function goes here
%   Detailed explanation goes here


str_folder="/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/dev/MultilevelEstimators/applications/SPDE/data/";
str_B="/scratch/leuven/330/vsc33032/Matrices/3d/";


level=maxlevel;
level_num=level;
HigherOrder=1;




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
         str_elements=str_folder+"Mesh/3d_beam/p_refinement/";
%          str_elements=str_folder+"GenerateMeshes_1_p_dev_beam/p_ref/Matlab/p_refinement/";

        
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





















file_nodes=strcat(str_elements,"Nodes_L",level,".txt");
files_elements=strcat(str_elements,"Elements_L",level,".txt");
Elements=dlmread(files_elements);
Nodes=dlmread(file_nodes);











switch level_num
    
    case 0
    increment=1;
    elementPlane='solid8';

    case 1
    increment=2;
    elementPlane='solid27_h';

    case 2
    increment=3;
    elementPlane='solid64_h';

    case 3
    increment=4;
    elementPlane='solid125_h';

    case 4
    increment=5;
    elementPlane='solid216_h';

end









%%%%%%%%%%

%% Parameters
Ek=ones(nElem,1).*30e6;         % Young's modulus in N/m^2
nu=0.15;        % Poisson ratio
el1=Elements(1,5:end);
nodesPos=Nodes(el1,:);
nod=nodesPos(:,2:end);
%     figure
%     for i=1:length(nod)
%         hold on
%         
%         plot3(nod(i,1),nod(i,2),nod(i,3),'*b')
%         text(nod(i,1),nod(i,2),nod(i,3),num2str((i)))
%         hold on        
%         
%         i=i+1;
%         
%     end



StepSize=0.1/10;


%% Create model
EltOpts.problem='3d';
EltOpts.nl='linear';
EltOpts.Coeffs=dlmread('coeffs_se125');
EltOpts.bendingmodes=false;
StrintPt='CC_SparseGrid_3d';
%StrintPt='gausslegendre_3d';

intPt=str2func( StrintPt );

EltOpts.nXi=8;
EltOpts.IntegrationFunc=intPt;
Types={1 elementPlane EltOpts};
Materials = cell(nElem,3);






%  for k =1:nElem, Materials(k,:)={k 'elastoplastic' {'isotropic' [E nu rho] 'dp' [alpha_dp k_dp(k,:)] 'multilinear' [1 k_dp(k,:)+E/300]}};end

 for k = 1:nElem, Materials(k,:) = {k 'linear' {'isotropic' [Ek(k,:) nu] }}; end

 
Sections=[1];
%




LeftNodes=selectnode(Nodes,-1e-6,-inf,-inf,1e-6,inf,inf);
RightNodes=selectnode(Nodes,max(Nodes(:,2))-1e-6,-inf,-inf,max(Nodes(:,2))+1e-6,inf,inf);
BottomNodes=selectnode(Nodes,-inf,-1e-6,-inf,inf,1e-6,inf);



% BottomNodes2=model.getNodes('lines',1);
% TopNodes=model.getNodes('lines',3);

LeftNodes=selectnode(Nodes,-1e-6,-inf,-inf,1e-6,inf,inf);
RightNodes=selectnode(Nodes,max(Nodes(:,2))-1e-6,-inf,-inf,max(Nodes(:,2))+1e-6,inf,inf);
DOF = getdof(Elements,Types);
sdof = [0.04;0.05;0.06;LeftNodes(:,1)+0.01;RightNodes(:,1)+0.01;LeftNodes(:,1)+0.02;RightNodes(:,1)+0.02;LeftNodes(:,1)+0.03;RightNodes(:,1)+0.03]; 
DOF = removedof(DOF,sdof);
sdog_algo=unique(floor(sdof));





%%%%%%%%%%%%%%%%%%%%%
% forceNodes=selectnode(Nodes,min(Nodes(:,2))-1e-6,min(Nodes(:,3))-1e-6,max(Nodes(:,4))-1e-6,max(Nodes(:,2))+1e-6,max(Nodes(:,3))+1e-6,max(Nodes(:,4))+1e-6);
% forceN=forceNodes(:,2:end);
% [srted,posid]=sortrows(forceN,[1 2]);
% n_nodes_x=length(uniquetol(Nodes(:,2),10^-7));
% n_nodes_y=length(uniquetol(Nodes(:,3),10^-7));
% n_nodes_z=length(uniquetol(Nodes(:,4),10^-7));
% forceN=forceNodes(posid,:);
% forceN=forceN(:,1);
% forceNodes=forceN+0.02;
%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
forceNodes=selectnode(Nodes,min(Nodes(:,2))-1e-6,max(Nodes(:,3))-1e-6,max(Nodes(:,4))./2-1e-6,max(Nodes(:,2))+1e-6,max(Nodes(:,3))+1e-6,max(Nodes(:,4))./2+1e-6);
forceN=forceNodes(:,2:end);
[srted,posid]=sortrows(forceN,[1 2]);
forceN=forceNodes(posid,:);
forceN=forceN(1:increment:end,:);
forceN=forceN(:,1);
forceNodes=forceN+0.02;
%%%%%%%%%%%%%%%%%

% forceNodes=selectnode(Nodes,max(Nodes(:,2))/2-1e-6,max(Nodes(:,3))./2-1e-6,max(Nodes(:,4))-1e-6,max(Nodes(:,2))/2+1e-6,max(Nodes(:,3))./2+1e-6,max(Nodes(:,4))+1e-6);
% forceN=forceNodes(:,1:end);
% forceN=forceN(:,1);
% forceNodes=forceN+0.02;



%% Load




%%%%%%%%%
%     Force=10000000*(0:1:1);
%     PLoad=Force;
%     
%  ar_x=[0:2*pi/(n_nodes_x-1):2*pi];
%  ar_y=[0:2*pi/(n_nodes_y-1):2*pi];
% vec_x=(-cos(ar_x)+1)/sum((-cos(ar_x)+1));
% vec_y=(-cos(ar_y)+1)/sum((-cos(ar_y)+1));
% mat=kron(vec_x,vec_y');
% vec=reshape(mat,[1,size(mat,2)*size(mat,1)]);
% P = nodalvalues(DOF,forceNodes,-vec);
%%%%%%%%%%%%%
    Force=10000000*(0:1:1);
    PLoad=Force;
ar=[0:2*pi/(length(forceNodes)-1):2*pi];
vec=(-cos(ar)+1)/sum((-cos(ar)+1));
P = -nodalvalues(DOF,forceNodes',vec);
%%%%%%%%%%
%     Force=10000000*(0:1:1);
%     PLoad=Force;
% P = -nodalvalues(DOF,forceNodes',-1);
%  figure;
%  plotnodes(Nodes);
% 
%  figure;
%  plotelem(Nodes,Elements,Types);

%% Compute equilibrium path
Options=struct;
Options.verbose=false;




[U,HistPar,critPnt,ElemCache,Divergence]=solver_nr(Nodes,Elements,Types,Sections,Materials,DOF,P,PLoad,Options,[],[],[],[],[]);


detjacEl=[];
for nel=1:nElem
    
  int=ElemCache{nel};  
 detjacEl(nel)=int.detJac;   
    
    
    
end

[val,pos]=uniquetol(detjacEl,10^-10);



if(length(pos)>1)
    fprintf("WARNING WARNING WARNING WARNINF\n") 
    fprintf("Multiple B matrices present\n")
end

dlmwrite(str_B+"detJac",val);

B=ElemCache{1}.B;
Bmat=B;

Order4From5=[1:7,9:11,13:15,17:19,21:29,37:43,45:47,49:51,53:55,57:65,73:79,81:83,85:87,89:91,93:101,109:115,117:119,121:123,125:127,129:137,181:187,189:191,193:195,197:199,201:209];
Order3from4=[1:6,8:9,11:12,14:15,17:20,26:31,33:34,36:37,39:40,42:45,51:56,58:59,61:62,64:65,67:70,101:106,108:109,111:112,114:115,117:120];
Order2from3=[1:5,7,9,11,14,17:21,23,25,27,29,49:53,55,57,59,61];
Order1from2=[1:4,19:22];
Collapse={};

Collapse{4}=Order4From5;
Collapse{3}=Order3from4;
Collapse{2}=Order2from3;
Collapse{1}=Order1from2;

        name="B_"+num2str((maxlevel+1));
       dlmwrite(str_B+name,Bmat);

      
for id=maxlevel:-1:1
    vec=Collapse{id};
    count=1;
    B_int=[];
        for len=1:length(vec)
       B_int(:,count*3-2:count*3)=Bmat(:,[vec(len)'.*3-2:vec(len)'.*3]);
       count=count+1;
        end
        name="B_"+num2str((id));
     Bmat=B_int;    

     dlmwrite(str_B+name,Bmat);

        
        
end


    
 QoINodes_3=selectnode(Nodes,-inf,min(Nodes(:,3))-1e-6,-1e-6+0.25,inf,min(Nodes(:,3))+1e-6,1e-6+0.25);



u_out_3=selectdof(DOF,QoINodes_3(:,1)+0.02)*U(:,end);



[s,id]=sort(QoINodes_3(:,2));
QoINodes_3=QoINodes_3(id,2);
u_out_3=u_out_3(id);   







end

