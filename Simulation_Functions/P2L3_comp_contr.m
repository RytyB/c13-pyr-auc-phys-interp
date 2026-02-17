function [pyr_fracs] = P2L3_comp_contr(parms, fdv)
    % Input: fit-parms and fdv for a 3PC model
    % Output: The relative contributions of each compartment to the pyruvate signal
    
    [EV, IV, vols] = P2L3(parms, fdv);
    
    pfb = sum( IV(1,:)*vols(1) );  % AUC Pyr IV
    pfe = sum( EV(1,:)*vols(2) );  % AUC Pyr EE
    pfc = sum( EV(3,:)*vols(3) );  % AUC Pyr C
    
    pyr_fracs = [pfb, pfe, pfc] / (pfb + pfe + pfc);
end

