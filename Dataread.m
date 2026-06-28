function [ThinElement, InitPosThin, ThickElement, InitPosThick] = Dataread(filename)
readfile = fopen(filename, 'r');
if readfile == -1
    error('Dataread:FileOpenFailed', 'Cannot open %s for reading.', filename);
end

cleanupObj = onCleanup(@() fclose(readfile));

skipLine(readfile);
counts = readNumericLine(readfile, 4);
NumThinEle = counts(1);
NumThinNode = counts(2);
NumThickEle = counts(3);
NumThickNode = counts(4);

skipLine(readfile);
ThinElement = zeros(NumThinEle, 2);
for i = 1:NumThinEle
    data = readNumericLine(readfile, 3);
    ThinElement(i, :) = data(2:3);
end

skipLine(readfile);
InitPosThin = zeros(NumThinNode, 1);
for i = 1:NumThinNode
    data = readNumericLine(readfile, 2);
    InitPosThin(data(1)) = data(2);
end

skipLine(readfile);
ThickElement = zeros(NumThickEle, 2);
for i = 1:NumThickEle
    data = readNumericLine(readfile, 3);
    ThickElement(i, :) = data(2:3);
end

skipLine(readfile);
InitPosThick = zeros(NumThickNode, 1);
for i = 1:NumThickNode
    data = readNumericLine(readfile, 2);
    InitPosThick(data(1)) = data(2);
end
end

function skipLine(fid)
line = fgetl(fid);
if ~ischar(line)
    error('Dataread:UnexpectedEOF', 'Unexpected end of geometry file.');
end
end

function data = readNumericLine(fid, minCount)
line = fgetl(fid);
if ~ischar(line)
    error('Dataread:UnexpectedEOF', 'Unexpected end of geometry file.');
end

data = sscanf(line, '%f').';
if numel(data) < minCount
    error('Dataread:BadLine', 'Bad numeric line in geometry file: %s', line);
end
end
