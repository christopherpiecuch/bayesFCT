%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function ClusMult = determine_clusters(COAST)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code last edited by CGP on 25 November 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ClusMult = determine_clusters(COAST)

% This script indicates which sea-level "cluster" a given tide-gauge site
% falls in (southeastern US or Caribbean+Central America+South America).

TARGET_COAST = [940 960];
BLOCK(find(ismember(COAST,TARGET_COAST))) = 1;
BLOCK(find(~ismember(COAST,TARGET_COAST))) = 2;

ClusMult = BLOCK'*BLOCK;
ClusMult(ClusMult==2)=0;
ClusMult(ClusMult~=0)=1;

return