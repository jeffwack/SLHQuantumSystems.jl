using SLHQuantumSystems

SQZ = squeezing_cavity(:A)
FCV = cavity(:B)

SYS = feedbackreduce(concatenate([SQZ,FCV],:sys),:A_Out,:B_In)

SS = StateSpace(SYS)
