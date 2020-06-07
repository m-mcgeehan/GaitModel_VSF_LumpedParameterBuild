plot(subj1.km.Rhip(:,2),'k'); hold on; plot(-1*subj1_calc.km.Rhip(:,2),'r'); 
plot(subj1.km.Lhip(:,2),'k'); hold on; plot(-1*subj1_calc.km.Lhip(:,2),'r'); 

plot(subj1.km.Rhip(:,1),'k'); hold on; plot(-1*subj1_calc.km.Rhip(:,1),'r'); 
plot(subj1.km.Lhip(:,1),'k'); hold on; plot(-1*subj1_calc.km.Lhip(:,1),'r'); 


plot(subj1.km.Rhip(:,3),'k'); hold on; plot(-1*subj1_calc.km.Rhip(:,3),'r'); 
plot(subj1.km.Lhip(:,3),'k'); hold on; plot(-1*subj1_calc.km.Lhip(:,3),'r'); 



plot(subj1.km.Rknee(:,2),'k'); hold on; plot(-1*subj1_calc.km.Rknee(:,2),'r'); 
plot(subj1.km.Lknee(:,2),'k'); hold on; plot(-1*subj1_calc.km.Lknee(:,2),'r'); 

a = [1 0 0;
    0 -1 0;
    0 0 1];
subj1.static.markers.pelvisOrigin = (.5*(subj1.static.markers.ASISmidpt+subj1.static.markers.PSISmidpt));
subj1.static.offsets.World2Pelvis = subj1.static.markers.pelvisOrigin;
subj1.static.offsets.World2Pelvis = subj1.static.offsets.World2Pelvis;
subj1.static.offsets.Pelvis2HJCR = subj1.static.markers.HJCR - subj1.static.offsets.World2Pelvis;
subj1.static.offsets.Pelvis2HJCL = a*subj1.static.offsets.Pelvis2HJCR;
subj1.static.offsets.HJCR2LegR =  subj1.static.markers.HJCR - subj1.static.v3dOffsets.LegRProxPos_GCS;
subj1.static.offsets.HJCL2LegL =  a*subj1.static.offsets.HJCR2LegR;
subj1.static.offsets.LegR2KJCR = [0 0 -subj1.static.segLength.legLengthR];
subj1.static.offsets.KJCR2ShankR = [0 1 -2.5];
subj1.static.offsets.LegL2KJCL = [0 0 -subj1.static.segLength.legLengthL];
subj1.static.offsets.KJCRLShankL = [0 -1 -2.5];
subj1.static.offsets.ShankL2AJCL = [-2 0 -subj1.static.segLength.shankLengthL];
shank2socket_offset = [-1.2 -1 -(subj1.static.segLength.shankLengthR+2.5)];
pylon_length = subj1.static.segLength.shankLengthL - (subj1.static.segLength.shankLengthR+2.5)+(subj1.static.segLength.footLZ - 8.7);



% subj1.static.offsets.LegR2KJCR = subj1_calcx.static.markers.KJCR - (subj1_calcx.static.markers.HJCR + subj1.static.offsets.HJCR2LegR);
% subj1.static.offsets.LegL2KJCL = subj1_calcx.static.markers.KJCL - (subj1_calcx.static.markers.HJCL +subj1.static.offsets.HJCL2LegL);
% subj1.static.offsets.KJCR2ShankR = subj1_calcx.static.markers.KJCR - subj1_calcx.static.v3dOffsets.ShankRProxPos_GCS;
% subj1.static.offsets.KJCL2ShankL = subj1_calcx.static.markers.KJCL - subj1_calcx.static.v3dOffsets.ShankLProxPos_GCS;
% subj1.static.offsets.ShankR2AJC_pros = subj1_calcx.static.markers.AJCR - subj1_calcx.static.v3dOffsets.ShankRProxPos_GCS;
% subj1.static.offsets.ShankL2AJCL = subj1_calcx.static.markers.AJCL - subj1_calcx.static.v3dOffsets.ShankLProxPos_GCS;
% 
% for 2segfoot
% subj1.static.markers.pelvisOrigin = (.5*(subj1.static.markers.ASISmidpt+subj1.static.markers.PSISmidpt));
% subj1.static.offsets.World2Pelvis = subj1.static.markers.pelvisOrigin;
% subj1.static.offsets.World2Pelvis = subj1.static.offsets.World2Pelvis;
% subj1.static.offsets.Pelvis2HJCR = (subj1.static.markers.HJCR - subj1.static.offsets.World2Pelvis);
% subj1.static.offsets.Pelvis2HJCL = (a*subj1.static.offsets.Pelvis2HJCR);
% subj1.static.offsets.HJCR2LegR =  (subj1.static.markers.HJCR - subj1.static.v3dOffsets.LegRProxPos_GCS);
% subj1.static.offsets.HJCL2LegL =  (a*subj1.static.offsets.HJCR2LegR);
% subj1.static.offsets.LegR2KJCR = [0 subj1.static.segLength.legLengthR 0];
% subj1.static.offsets.KJCR2ShankR = [0 0 2.5];
% subj1.static.offsets.LegL2KJCL = [0 subj1.static.segLength.legLengthL 0];
% subj1.static.offsets.KJCRLShankL = [0 0 2.5];
% subj1.static.offsets.ShankL2AJCL = [0 0 subj1.static.segLength.shankLengthL];
% shank2socket_offset = [1.2 1 (subj1.static.segLength.shankLengthR+2.5)];
% pylon_length = subj1.static.segLength.shankLengthL - (subj1.static.segLength.shankLengthR+2.5)+(subj1.static.segLength.footLZ - 8.7);
% footL_anatomic_offset = [1.0 0 -1.5];

