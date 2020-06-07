addpath(genpath('matlab scripts'),...
        genpath('stlTools'),...
        genpath('Subject Bone Geometry')); %Add dependencies
oldpath = cd;
cd('Subject Bone Geometry');

disp('Scaling bone geometry and assembling model...');
%% Read in generic geometries
[v, f, n, name] = stlRead('FemurR.stl');
FemurR_sf = dynamic.static.segLength.legLengthR*100/41.58; %Create uniform scaling factor between the generic bone length and the bone length from MoCap
fv.vertices = (v*FemurR_sf)/100; %Convert units
fv.faces = f;
dynamic.boneGeom.femurR = fv; 
stlWrite('FemurR_scaled.stl',fv); 
clear v f n name fv

[v, f, n, name] = stlRead('FemurL.stl');
FemurL_sf = dynamic.static.segLength.legLengthL*100/41.58;
fv.vertices = (v*FemurL_sf)/100;
fv.faces = f;
dynamic.boneGeom.femurL = fv;
stlWrite('FemurL_scaled.stl',fv); 
clear v f n name fv

[v, f, n, name] = stlRead('ShankL.stl');
ShankL_sf = dynamic.static.segLength.shankLengthL*100/37.25;
fv.vertices = (v*ShankL_sf)/100;
fv.faces = f;
dynamic.boneGeom.shankL = fv;
stlWrite('ShankL_scaled.stl',fv); 
clear v f n name fv

[v, f, n, name] = stlRead('ShankR_cut.stl');
ShankR_cut_sf = (.66*(dynamic.static.segLength.shankLengthL*100))/19.0;
fv.vertices = [v(:,1)*ShankL_sf v(:,2)*ShankL_sf v(:,3)*ShankR_cut_sf]/100;
fv.faces = f;
dynamic.boneGeom.shankR_cut = fv;
stlWrite('ShankR_cut_scaled.stl',fv); 
clear v f n name fv

[v, f, n, name] = stlRead('Pelvis.stl'); 
% pelvis_sf = (dynamic.static.segLength.ASISYdist/25.7+ dynamic.static.segLength.PSISYdist/8.63)/2;
pelvis_sfX = dynamic.static.segLength.pelvisX/14.5;
pelvis_sfY = (dynamic.static.segLength.ASISYdist/25.8 + dynamic.static.segLength.PSISYdist/9.14)/2;
pelvis_sfZ = (pelvis_sfX+pelvis_sfY)/2;
fv.vertices = [v(:,1)*pelvis_sfX v(:,2)*pelvis_sfY v(:,3)*pelvis_sfZ]/100;
fv.faces = f;
dynamic.boneGeom.pelvis = fv;
stlWrite('Pelvis_scaled.stl',fv); 
clear v f n name fv

[v, f, n, name] = stlRead('HAT.stl');
HAT_sf = (pelvis_sfX + pelvis_sfY + pelvis_sfZ)/3;
fv.vertices = (v*HAT_sf)/100;
fv.faces = f;
dynamic.boneGeom.HAT = fv;
stlWrite('HAT_scaled.stl',fv); 
clear v f n name fv


[v, f, n, name] = stlRead('FootL.stl');
FootL_sf = dynamic.static.segLength.footLX/21.79;
fv.vertices = (v*FootL_sf)/100;
fv.faces = f;
stlWrite('FootL_scaled.stl',fv); 
clear v f n name fv

disp('Scaling bone geometry and assembling model...Done');

cd(oldpath);
