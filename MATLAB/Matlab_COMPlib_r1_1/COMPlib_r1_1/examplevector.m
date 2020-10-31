%
%  Matlab script file: examplevector.m
%
%  It provides several string vectors containing the names of certain
%  groups of the benchmark examples in COMPleib version 1.1.
%
%  Note, this script file can be used in combination with 
%  the matlab function COMPleib.m, i.e. 
%
%  >> examplevector;   %--- script file defining ex.-vectors of COMPlib
%  >> Set_of_Ex=Ex;    %--- string vector of all examples of COMPleib 1.0
%  >>                  %    (see examplevector.m for def. of the string vector Ex)
%  >> [No_Ex,dummy]=size(Set_of_Ex);
%  >> for i=1:No_Ex
%  >>     %--- load COMPleib test examples
%  >>     [A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny]=COMPleib(Set_of_Ex(i,:));
%  >>      ......  (individual matlab code, i.e. define a NSDP from COMPleib
%  >>      ......   and solve it by your NSDP-solver)
%  >> end
%
%  For more details, see [1], [2], [3].
%
%     References: 
%     [1] COMPleib: COnstrained Matrix-optimization Problem library --
%         a collection of test examples for nonlinear semidefinite
%         programs, control system design and related problems
%         F. Leibfritz, Tech.-Report, 2003
%
%     [2] COMPleib 1.0 -- User manual and quick reference
%         F. Leibfritz and W. Lipinski, Tech.-Report, 2004
%
%     [3] Description of the benchmark examples in COMPl{_e}ib 1.0
%         F. Leibfritz and W. Lipinski, Tech.-Report, 2003
% 
%
%  Author: Friedemann Leibfritz and Wojtek Lipinski
%
%    Note: Wojtek is a diploma student at the University of Trier.
%          I want to thank Wojtek since he has done a very good job
%          in the COMPleib project.
%
%    Date: 23.12.2003
%          09.03.2005 (Leibfritz)

ExAC  =strvcat('AC1' ,'AC2' ,'AC3' ,'AC4' ,'AC5' ,'AC6' ,'AC7' ,'AC8' ,'AC9' ,'AC10','AC11','AC12','AC13','AC14','AC15','AC16','AC17','AC18'); 
ExHE  =strvcat('HE1' ,'HE2' ,'HE3' ,'HE4' ,'HE5' ,'HE6' ,'HE7'); 
ExJE  =strvcat('JE1' ,'JE2' ,'JE3');     
ExREA =strvcat('REA1','REA2','REA3','REA4');
ExDIS =strvcat('DIS1','DIS2','DIS3','DIS4','DIS5'); 
ExWEC =strvcat('WEC1','WEC2','WEC3');
ExEB  =strvcat('EB1' ,'EB2' ,'EB3' ,'EB4' ,'EB5' ,'EB6' );
ExTF  =strvcat('TF1' ,'TF2' ,'TF3' ); 
ExNN  =strvcat('NN1' ,'NN2' ,'NN3' ,'NN4' ,'NN5' ,'NN6' ,'NN7' ,'NN8' ,'NN9' ,'NN10','NN11','NN12','NN13','NN14','NN15','NN16','NN17','NN18');
% other:
Exsof =strvcat('TG1' ,'AGS' ,'HF1' ,'BDT1','BDT2','MFP' ,'UWV' ,'IH'  ,'CSE1','CSE2','PAS' ,'PSM' ,'TL'  ,'CDP' );

% --> all pure SOF examples:
ExSOF =strvcat(ExAC,ExHE,ExJE,ExREA,ExDIS,'TG1','AGS',ExWEC,'HF1','BDT1','BDT2','MFP','UWV','IH','CSE1','CSE2',ExEB,'PAS',ExTF,'PSM','TL','CDP',ExNN); 

%-------------------------------------------------------------------------------
% heat flow 2D models (small POD-approx., dense):
ExHF2D_cd_pod=strvcat('HF2D_CD4','HF2D_CD5','HF2D_CD6');
ExHF2D_is_pod=strvcat('HF2D_IS5','HF2D_IS6','HF2D_IS7','HF2D_IS8');
ExHF2D_pod=strvcat(ExHF2D_is_pod,ExHF2D_cd_pod,'HF2D10','HF2D11','HF2D12','HF2D13','HF2D14','HF2D15','HF2D16','HF2D17','HF2D18');

%-------------------------------------------------------------------------------
% heat flow 2D models (medium scale, sparse):
ExHF2D_is_medium=strvcat('HF2D_IS1_M361','HF2D_IS1_M529','HF2D_IS2_M361','HF2D_IS2_M529','HF2D_IS3_M256','HF2D_IS3_M484','HF2D_IS4_M256','HF2D_IS4_M484');
ExHF2D_cd_medium=strvcat('HF2D_CD1_M256','HF2D_CD1_M484','HF2D_CD2_M256','HF2D_CD2_M484','HF2D_CD3_M324','HF2D_CD3_M576');
ExHF2D_medium=strvcat(ExHF2D_is_medium,ExHF2D_cd_medium,'HF2D1_M316','HF2D1_M541','HF2D2_M316','HF2D2_M541','HF2D5_M289','HF2D5_M529','HF2D6_M289','HF2D6_M529','HF2D9_M256','HF2D9_M484');

%-------------------------------------------------------------------------------
% heat flow 2D models (large scale, sparse):
ExHF2D_large=strvcat('HF2D_IS1','HF2D_IS2','HF2D_IS3','HF2D_IS4','HF2D_CD1','HF2D_CD2','HF2D_CD3','HF2D1','HF2D2','HF2D3','HF2D4','HF2D5','HF2D6','HF2D7','HF2D8','HF2D9');


%-------------------------------------------------------------------------------
% all heat flow 2D models:
ExHF2D=strvcat(ExHF2D_large,ExHF2D_medium,ExHF2D_pod);


%-------------------------------------------------------------------------------
ExCM  =strvcat('CM1' ,'CM2' ,'CM3' ,'CM4' ,'CM5' ,'CM6' ); 
ExCM_IS  =strvcat('CM1_IS' ,'CM2_IS' ,'CM3_IS' ,'CM4_IS' ,'CM5_IS' ,'CM6_IS' );
% other:
Ex2om =strvcat('TMD' ,'FS'  ,'DLR1','DLR2','DLR3','ISS1','ISS2','CBM' ,'LAH' );

% --> all second order models: 
Ex2OM =strvcat(ExCM,ExCM_IS,Ex2om);

%-------------------------------------------------------------------------------
% reduced order control systems:
ExROC =strvcat('ROC1','ROC2','ROC3','ROC4','ROC5','ROC6','ROC7','ROC8','ROC9','ROC10'); 

%===============================================================================
% ==> all examples of COMPleib 1.1:
Ex=strvcat(ExSOF,ExHF2D,Ex2OM,ExROC);

% all examples without "Exsparse_large":
Ex_no_sl=strvcat(ExSOF(1:80,:),ExHF2D_medium,ExHF2D_pod,Ex2OM,ExROC);

%===============================================================================
% all examples of order order >= 50:
Exlarge=strvcat('AC10','HF1','BDT2','CSE2','EB6','TL','CDP','NN18',ExHF2D_large,ExHF2D_medium,ExCM(2:6,:),ExCM_IS(2:6,:),'ISS1','ISS2','CBM');
% all examples of order < 50:
Exsmall=strvcat(ExAC(1:9,:),ExAC(11:18,:),ExHE,ExJE,ExREA,ExDIS,'TG1','AGS',ExWEC,'BDT1','MFP','UWV','IH','CSE1',ExEB(1:5,:),'PAS',ExTF,'PSM',ExNN(1:17,:),ExHF2D_pod,'CM1','CM1_IS','TMD','FS','DLR1','DLR2','DLR3','LAH',ExROC);

% all dense examples:
Exdense  =strvcat(ExAC(1:9,:),'AC11','AC12',ExAC(15:18,:),ExHE,'JE2','JE3',ExREA,ExDIS,'TG1',ExWEC,'MFP','UWV','PAS',ExTF,'PSM','TL',ExNN(1:17,:),ExHF2D_pod,'TMD','FS','DLR1',ExROC);     
% all partly sparse examples:
Expsparse=strvcat('JE1',ExCM,ExCM_IS,'CBM','LAH');
% all sparse examples, but without the ones in "Exsparse_large":
Exsparse =strvcat(ExHF2D_medium,'AC10','AC13','AC14','AGS','HF1','BDT1','BDT2','IH','CSE1','CSE2',ExEB,'CDP','DLR2','DLR3','ISS1','ISS2');
% sparse examples of very high order saved in a sparse format:
Exsparse_large =strvcat('NN18',ExHF2D_large);
