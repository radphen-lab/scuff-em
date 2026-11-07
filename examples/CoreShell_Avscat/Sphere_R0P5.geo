SetFactory("OpenCASCADE");
Sphere(1) = {0,0,0,0.5};
Mesh.CharacteristicLengthMin = 0.11;
Mesh.CharacteristicLengthMax = 0.11;
Mesh.MshFileVersion = 2.2;
Mesh.ElementOrder = 1;
Physical Surface(1) = {1};
