function [MotorSK, ForceMotorVector] = MotorSKDD(MotorBindSite, MotorWorkState, SwingDis, initPosition)
global NumThinNode NumThickNode MotorStiff

nodeNum = NumThinNode + NumThickNode;
MotorSK = zeros(nodeNum, nodeNum);
ForceMotorVector = zeros(nodeNum, 1);

for i = 1:NumThickNode
    if abs(MotorWorkState(i)) == 2
        thinNode = MotorBindSite(i);
        thickNode = NumThinNode + i;
        motorExtension = SwingDis(i) + initPosition(thinNode) - initPosition(thickNode);

        MotorSK(thinNode, thinNode) = MotorSK(thinNode, thinNode) + MotorStiff;
        MotorSK(thinNode, thickNode) = MotorSK(thinNode, thickNode) - MotorStiff;
        MotorSK(thickNode, thinNode) = MotorSK(thickNode, thinNode) - MotorStiff;
        MotorSK(thickNode, thickNode) = MotorSK(thickNode, thickNode) + MotorStiff;

        ForceMotorVector(thinNode) = ForceMotorVector(thinNode) + MotorStiff * motorExtension;
        ForceMotorVector(thickNode) = ForceMotorVector(thickNode) - MotorStiff * motorExtension;
    end
end
end

