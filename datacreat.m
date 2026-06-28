function datacreat(filename)
if nargin < 1 || isempty(filename)
    filename = fullfile(fileparts(mfilename('fullpath')), 'Data4readRigid.txt');
end

NumThinEle = 222;
NumThinNode = 223;
NumThickEle = 76;
NumThickNode = 77;

PosThinNode = (0:NumThinNode - 1) * 5.5;
PosThickNode = (NumThickNode - 1:-1:0) * 14.3;

Inputfile = fopen(filename, 'w');
if Inputfile == -1
    error('datacreat:FileOpenFailed', 'Cannot open %s for writing.', filename);
end

fprintf(Inputfile, 'Numbers of elements or nodes(thin element, thin nodes, thick elements, thick nodes)\n');
fprintf(Inputfile, '%8d %8d %8d %8d\n', NumThinEle, NumThinNode, NumThickEle, NumThickNode);
fprintf(Inputfile, 'Information of elements which belong to thin filament\n');
ThinEleID = 1:NumThinEle;
fprintf(Inputfile, '%8d %8d %8d\n', [ThinEleID; ThinEleID; ThinEleID + 1]);
fprintf(Inputfile, 'Information of nodes on thin filament\n');
fprintf(Inputfile, '%8d %10.5f\n', [1:NumThinNode; PosThinNode]);
fprintf(Inputfile, 'Information of elements which belong to thick filament\n');
ThickEleID = 1:NumThickEle;
fprintf(Inputfile, '%8d %8d %8d\n', [ThickEleID; ThickEleID; ThickEleID + 1]);
fprintf(Inputfile, 'Information of nodes on thick filament\n');
fprintf(Inputfile, '%8d %10.5f\n', [1:NumThickNode; PosThickNode]);
fclose(Inputfile);
end
