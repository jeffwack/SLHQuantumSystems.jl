using SLHQuantumSystems

SQZ = squeezing_cavity("A")
FCV = cavity("B")

SYS = concatenate([SQZ,FCV],"sys")
SYS = feedbackreduce(SYS,"A_out","B_in")

SS = StateSpace(SYS)
