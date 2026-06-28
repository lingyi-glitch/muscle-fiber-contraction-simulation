function InitSK = InitSKDD(ThinElement, ThickElement)
global NumThinEle NumThinNode NumThickEle NumThickNode StiffThin StiffThick kthick

nodeNum = NumThinNode + NumThickNode;
InitSK = zeros(nodeNum, nodeNum);

for i = 1:NumThinEle
    n1 = ThinElement(i, 1);
    n2 = ThinElement(i, 2);
    InitSK(n1, n1) = InitSK(n1, n1) + StiffThin;
    InitSK(n1, n2) = InitSK(n1, n2) - StiffThin;
    InitSK(n2, n1) = InitSK(n2, n1) - StiffThin;
    InitSK(n2, n2) = InitSK(n2, n2) + StiffThin;
end

for i = 1:NumThickEle
    n1 = ThickElement(i, 1) + NumThinNode;
    n2 = ThickElement(i, 2) + NumThinNode;
    InitSK(n1, n1) = InitSK(n1, n1) + StiffThick;
    InitSK(n1, n2) = InitSK(n1, n2) - StiffThick;
    InitSK(n2, n1) = InitSK(n2, n1) - StiffThick;
    InitSK(n2, n2) = InitSK(n2, n2) + StiffThick;
end

fixedNode = NumThinNode + NumThickNode;
InitSK(fixedNode, fixedNode) = InitSK(fixedNode, fixedNode) + kthick;
end

